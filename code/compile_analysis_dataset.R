#!/usr/bin/env Rscript

# here are our standard example cases:
#   single antibody:  VRC07-523-LS
#   double cocktail:  VRC07-523LS + PGT121.BIJ4141LS
#   triple cocktail:  VRC07-523LS + PGT121.BIJ414LS + PGDM1400LS


# ---------------------------------------------------------------------------- #
# STEP -1:  prepare our environment
# ---------------------------------------------------------------------------- #
library(seqinr)
path.home <- "/home"

# antibody names are passed to docker container at run time as
# environment variable Nab, which is a semicolon-separated list
antibody_string <- Sys.getenv("nab")
sensitivity.threshold <- Sys.getenv("sens_thresh")
multiple.sensitivity.threshold <- Sys.getenv("multsens_nab")
antibodies <- strsplit(antibody_string, split = ";")[[1]]
filename <- Sys.getenv("nab_str")

# ---------------------------------------------------------------------------- #
# STEP 0:  load and prepare our data
# ---------------------------------------------------------------------------- #


# define our working directories
path.lib <- file.path(path.home, "lib")
path.data <- file.path(path.home, "dat")
path.data.catnap <- file.path(path.data, "catnap")
path.data.analysis <- file.path(path.data, "analysis")
path.out <- file.path(path.home, "output")

source(file.path(path.lib, "utils.R"))
opts <- get_global_options()
# load data
data.assay <- read.table(file.path(path.data.catnap, "assay.txt"), header=T, sep="\t", quote="\"")
data.viruses <- read.table(file.path(path.data.catnap, "viruses.txt"), header=T, sep="\t", quote="\"")
data.abs <- read.table(file.path(path.data.catnap, "abs.txt"), header=T, sep="\t", quote="\"")

# source our function library
source(file.path(path.lib, "multi_ab.Rlib"))

# load and process virus info and sequences
data.seq <- read.fasta(file.path(path.data.catnap, "virseqs_aa.fasta"), seqtype="AA")
seqname.full <- names(data.seq)
header.info <- strsplit(names(data.seq), split=".", fixed=T)
subtype <- unlist(lapply(header.info, function(x) return(x[1])))
country <- unlist(lapply(header.info, function(x) return(x[2])))
year <- unlist(lapply(header.info, function(x) return(x[3])))
seqname.db <- unlist(lapply(header.info, function(x) return(x[4])))

# let's filter out outlier assay results(with IC50 of ">1")
data.assay.reduced <- data.assay[data.assay$IC50 != ">1", ]

# determine which sequences were assayed for all antibodies in the cocktail
ab.seqs <- list()
for(ab.index in 1:length(antibodies)) {
  ab.seqs[[ab.index]] <- unique(as.character(data.assay.reduced[data.assay.reduced[, 1] == antibodies[ab.index], 2]))
}
seqs.include <- Reduce(intersect, ab.seqs)

# let's define the data relevant to the selected antibod(y|ies)
seqs.selected <- data.seq[seqname.full[seqname.db %in% seqs.include]]
seqname.selected.full <- names(seqs.selected)
seqname.selected.db <- seqname.db[seqname.full %in% seqname.selected.full]
subtype.selected <- subtype[seqname.full %in% seqname.selected.full]
country.selected <- country[seqname.full %in% seqname.selected.full]


# ---------------------------------------------------------------------------- #
# STEP 1:  process our IC50/IC80 readouts
# ---------------------------------------------------------------------------- #


# initialize our readout vectors
imputed.ic50 <- list()
imputed.ic80 <- list()
censored.ic50 <- list()
censored.ic80 <- list()
for(ab.tmp in antibodies) {
    ab_tmp_str <- paste0("nab_", ab.tmp)
  imputed.ic50[[ab_tmp_str]] <- rep(NA, length(seqs.selected))
  imputed.ic80[[ab_tmp_str]] <- rep(NA, length(seqs.selected))
  censored.ic50[[ab_tmp_str]] <- rep(0, length(seqs.selected))
  censored.ic80[[ab_tmp_str]] <- rep(0, length(seqs.selected))
}

# collect and process our imputed/censored information
for(ab.tmp in antibodies) {
    ab_tmp_str <- paste0("nab_", ab.tmp)
  for(seqname.index in 1:length(seqname.selected.db)) {

    # isolate our readouts of interest
    seqname.tmp <- seqname.selected.db[seqname.index]
    data.assay.reduced.tmp <- data.assay.reduced[data.assay.reduced$Antibody == ab.tmp & data.assay.reduced$Virus == seqname.tmp, 5:6]

    # IC50:  let's confirm that we actually have readout data for this sequence
    if(length(data.assay.reduced.tmp[data.assay.reduced.tmp[, 1] != "", 1]) > 0) {

      # make the binary notation that we have a right-censored variable
      if(">" %in% substr(data.assay.reduced.tmp[, 1], 1, 1)) {
        censored.ic50[[ab_tmp_str]][seqname.index] <- 1
      }

      # merge all of our readouts, both censored and numeric
      imputed.ic50[[ab_tmp_str]][seqname.index] <- merge.readouts(data.assay.reduced.tmp[, 1])

    # if no assay readouts exist, "Mark it zero"
    } else {
      censored.ic50[[ab_tmp_str]][seqname.index] <- NA
      imputed.ic50[[ab_tmp_str]][seqname.index] <- NA
    }

    # IC80:  let's confirm that we actually have readout data for this sequence
    if(length(data.assay.reduced.tmp[data.assay.reduced.tmp[, 2] != "", 2])) {

      # make the binary notation that we have a right-censored variable
      if(">" %in% substr(data.assay.reduced.tmp[, 2], 1, 1)) {
        censored.ic80[[ab_tmp_str]][seqname.index] <- 1
      }

      # merge all of our readouts, both censored and numeric
      imputed.ic80[[ab_tmp_str]][seqname.index] <- merge.readouts(data.assay.reduced.tmp[, 2])

    # if no assay readouts exist, "Mark it zero"
    } else {
      censored.ic80[[ab_tmp_str]][seqname.index] <- NA
      imputed.ic80[[ab_tmp_str]][seqname.index] <- NA
    }
  }
}

# combine results into new dataframe
readouts <- data.frame(seq.id.catnap=seqname.selected.db)
for(ab.tmp in antibodies) {
    ab_tmp_str <- paste0("nab_", ab.tmp)
  readouts.tmp <- data.frame(censored.ic50[[ab_tmp_str]], censored.ic80[[ab_tmp_str]],
                              imputed.ic50[[ab_tmp_str]], imputed.ic80[[ab_tmp_str]])
  names(readouts.tmp) <- c(paste0(ab_tmp_str, ".ic50.censored"),
                             paste0(ab_tmp_str, ".ic80.censored"),
                             paste0(ab_tmp_str, ".ic50.imputed"),
                             paste0(ab_tmp_str, ".ic80.imputed"))
  readouts <- data.frame(readouts, readouts.tmp)
}

# if we are prespecifying more than one antibody, then let's predict the
# combined outcomes
if(length(antibodies) > 1) {
  # use additive method to determine predicted combinations of IC50/IC80
  #(i.e., the "Quantitative 1" endpoint)
  readouts$pc.ic50 <- apply(readouts[, grep("ic50.imputed", names(readouts), fixed=T)], 1, wagh.additive.method)
  readouts$pc.ic80 <- apply(readouts[, grep("ic80.imputed", names(readouts), fixed=T)], 1, wagh.additive.method)
} else { # if only one antibody, define them as the single-ab version
    readouts$pc.ic50 <- readouts[, grep("ic50.imputed", names(readouts), fixed = TRUE)]
    readouts$log10.pc.ic50 <- log10(readouts$pc.ic50)
    readouts$pc.ic80 <- readouts[, grep("ic80.imputed", names(readouts), fixed = TRUE)]
    readouts$log10.pc.ic80 <- log10(readouts$pc.ic80)
}
readouts$log10.pc.ic50 <- log10(readouts$pc.ic50)
readouts$log10.pc.ic80 <- log10(readouts$pc.ic80)
# calculate IIP(a.k.a. the "Quantitive 2" endpoint)
iip.c <- 10
iip.m <- log10(4) /(readouts$log10.pc.ic80 - readouts$log10.pc.ic50)
iip.f.c <-(iip.c ^ iip.m) /((readouts$pc.ic50 ^ iip.m) +(iip.c ^ iip.m))
iip.f.c[iip.f.c >= 1] <- 1 - .Machine$double.neg.eps
readouts$iip <- (-1) * log10(1 - iip.f.c)

# derive the "Dichotomous 1" endpoint(i.e., is the PC IC50 less than the
# sensitivity cutoff?)
readouts$dichotomous.1 <- as.numeric(readouts$pc.ic50 < sensitivity.threshold)

# derive the "Dichotomous 2" endpoint(i.e., is the imputed IC50 less than
# the sensitivity threshold for at least a user-select number of antibodies
min.sensitive.abs <- min(c(multiple.sensitivity.threshold, length(antibodies)))
readouts$dichotomous.2 <- as.numeric(apply(readouts[, grep("ic50.imputed", names(readouts), fixed=TRUE), drop = FALSE] < sensitivity.threshold, 1, sum) >= min.sensitive.abs)


# ---------------------------------------------------------------------------- #
# STEP 2:  take care of business
# ---------------------------------------------------------------------------- #


# let's start off by making our HXB2 map
hxb2.seq <- toupper(paste(data.seq$B.FR.1983.HXB2.K03455, collapse=""))
hxb2.map <- mk.hxb2.map(hxb2.seq)

# generate our AA residue information(and remove invariant columns)
env.aa.indicators <- as.data.frame(aa.char.2.aa.binary(seqs.selected, seqname.selected.full, hxb2.map))
env.aa.indicators.reduced <- env.aa.indicators[, apply(env.aa.indicators, 2, var) != 0]

# create vectors of PNGS indicator values
env.sequon.indicators <- as.data.frame(aa.char.2.sequon.indicators(seqs.selected, seqname.selected.full, hxb2.map))
env.sequon.indicators.reduced <- env.sequon.indicators[, colSums(env.sequon.indicators) > 0]

# determine our viral geometry measures(and remove invariant columns)
env.virgeom <- aa.char.2.vir.geom(seqs.selected, seqname.selected.full, hxb2.map)
env.virgeom.reduced <- env.virgeom[, apply(env.virgeom, 2, var) != 0]

# assemble our data set
data.final <- data.frame(seq.id.lanl=seqname.selected.full,
                          seq.id.catnap=seqname.selected.db,
                          readouts[, -1],
                          as.data.frame(country.2.geo.reg(country.selected)),
                          as.data.frame(bin.subtype(subtype.selected)),
                          env.aa.indicators.reduced,
                          env.sequon.indicators.reduced,
                          env.virgeom.reduced)

# remove features pertaining to insertions in HXB2
filter.insertions <- rep (TRUE, ncol (data.final))
for (var.index in 1:ncol(data.final)) {
  varname.components <- unlist (strsplit (names (data.final)[var.index], split=".", fixed=T))
  if (varname.components[1] == "hxb2" & grepl ('[a-z]', varname.components[2])) {
    filter.insertions[var.index] <- FALSE
  }
}
data.final <- data.final[ , filter.insertions]

# name our outfile and save
current_date <- as.Date(Sys.getenv('current_date'), "%d%b%Y")
# filename <- paste0("slapnap_", paste(gsub("/", "-", antibodies), collapse="_"), "_", format(current_date, "%d%b%Y"), ".csv")
final_filename <- paste0("slapnap_", filename, "_", format(current_date, "%d%b%Y"), ".csv")
setwd(path.data.analysis)
write.csv(data.final, file=final_filename, row.names=F)

# save missing data stats for report compilation later
nprevious <- length(data.final[,1])
saveRDS(nprevious, "/home/slfits/nprevious.rds")

# first covariate column
min_cov_col_idx <- min(grep("geographic", colnames(data.final)))
ncol_data_final <- ncol(data.final)

# number with complete sequence/geog information
complete_features_idx <- complete.cases(data.final[ , min_cov_col_idx:ncol_data_final])
ncomplete_features <- sum(complete_features_idx)

# number with complete IC50 + complete features
ncomplete_ic50 <- sum(!is.na(data.final$pc.ic50[complete_features_idx]))

# number with complete IC80 + complete features
ncomplete_ic80 <- sum(!is.na(data.final$pc.ic80[complete_features_idx]))

# number with complete IC50 and IC80 + complete features
ncomplete_ic5080 <- sum(!is.na(data.final$pc.ic50[complete_features_idx]) & !is.na(data.final$pc.ic80[complete_features_idx]))

saveRDS(ncomplete_features, "/home/slfits/ncomplete_features.rds")
saveRDS(ncomplete_ic50, "/home/slfits/ncomplete_ic50.rds")
saveRDS(ncomplete_ic80, "/home/slfits/ncomplete_ic80.rds")
saveRDS(ncomplete_ic5080, "/home/slfits/ncomplete_ic5080.rds")
# ---------------------------------------------------------------------------- #
#                                    - 30 -
# ---------------------------------------------------------------------------- #
