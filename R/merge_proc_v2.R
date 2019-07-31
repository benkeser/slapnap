#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# source ("/home/cmagaret/compbio/projectus/cam117_multiAb/bin/02_procdata/merge_proc_v2.R")
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# here are our standard example cases:
#   double cocktail:  VRC07-523LS + PGT12.BIJ4141LS 
#   triple cocktail:  VRC07-523LS + PGT121.BIJ414LS + PGDM1400LS


# ---------------------------------------------------------------------------- #
# STEP -1:  prepare our environment
# ---------------------------------------------------------------------------- #


# refresh our workspace
rm (list = ls ())

# what libraries do we need to survive?
library (seqinr)
#library (gtools)

# CONFIG:  define our home directory
path.home <- "/home"

# CONFIG:  which antibodies are we interested in?
#antibodies <- c ("VRC07-523-LS", "PGT121")
# antibodies <- c ("VRC07-523-LS", "PGT121", "PGDM1400")

# antibody names are passed to docker container at run time as 
# environment variable Nab, which is a semicolon-separated list
antibody_string <- Sys.getenv("Nab")
antibodies <- strsplit(antibody_string, split = ";")[[1]]

# ---------------------------------------------------------------------------- #
# STEP 0:  load and prepare our data
# ---------------------------------------------------------------------------- #


# define our working directories
path.lib <- file.path (path.home, "lib")
path.data <- file.path (path.home, "dat")
path.data.catnap <- file.path (path.data, "catnap")
path.data.analysis <- file.path (path.data, "analysis")

# load data
data.assay <- read.table (file.path (path.data.catnap, "assay.txt"), header = TRUE, sep="\t", quote="\"")
data.viruses <- read.table (file.path (path.data.catnap, "viruses.txt"), header = TRUE, sep="\t", quote="\"")
data.abs <- read.table (file.path (path.data.catnap, "abs.txt"), header = TRUE, sep="\t", quote="\"")

# source our function library
source (file.path (path.lib, "multi_ab_v2.Rlib"))

# load and process virus info and sequences
data.seq <- read.fasta (file.path (path.data.catnap, "virseqs_aa.fasta"), seqtype="AA")
seqname.full <- names (data.seq)
header.info <- strsplit (names (data.seq), split=".", fixed = TRUE)
subtype <- unlist (lapply (header.info, function (x) return (x[1])))
country <- unlist (lapply (header.info, function (x) return (x[2])))
seqname.db <- unlist (lapply (header.info, function (x) return (x[4]))) # entry 3 appears now to be a year

# prep information for our antibody combinations
#ab.combinations <- apply (permutations (n=length (antibodies), r=length (antibodies), v=antibodies), 1, paste0, collapse="+")

# determine which sequences were assayed for all antibodies in the cocktail
ab.seqs <- list ()
for (ab.index in 1:length (antibodies)) {
  ab.seqs[[ab.index]] <- unique (as.character (data.assay[data.assay[, 1] == antibodies[ab.index], 2]))
}
seqs.include <- Reduce (intersect, ab.seqs)

# let's define the data relevant to the selected antibod(y|ies)
seqs.selected <- data.seq[seqname.full[seqname.db %in% seqs.include]]
seqname.selected.full <- names (seqs.selected)
seqname.selected.db <- seqname.db[seqname.full %in% seqname.selected.full]
subtype.selected <- subtype[seqname.full %in% seqname.selected.full]
country.selected <- country[seqname.full %in% seqname.selected.full]


# ---------------------------------------------------------------------------- #
# STEP 1:  process our IC50/IC80 readouts
# ---------------------------------------------------------------------------- #


# initialize our readout vectors
imputed.ic50 <- list ()
imputed.ic80 <- list ()
censored.ic50 <- list ()
censored.ic80 <- list ()
for (ab.tmp in antibodies) {
  imputed.ic50[[ab.tmp]] <- rep (NA, length (seqs.selected))
  imputed.ic80[[ab.tmp]] <- rep (NA, length (seqs.selected))
  censored.ic50[[ab.tmp]] <- rep (0, length (seqs.selected))
  censored.ic80[[ab.tmp]] <- rep (0, length (seqs.selected))
}

# collect our imputed/censored information
for (ab.tmp in antibodies) {
  for (seqname.index in 1:length (seqname.selected.db)) {

    # isolate our readouts of interest
    seqname.tmp <- seqname.selected.db[seqname.index]
    data.assay.tmp <- data.assay[data.assay$Antibody == ab.tmp & data.assay$Virus == seqname.tmp, 5:6]

    # we need to make sure that we actually have readout data for this sequence
    if (length (data.assay.tmp[data.assay.tmp[, 1] != "", 1])) {
      if (">" %in% substr (data.assay.tmp[, 1], 1, 1)) {
        censored.ic50[[ab.tmp]][seqname.index] <- 1
      }
      imputed.ic50[[ab.tmp]][seqname.index] <- merge.readouts (data.assay.tmp[, 1])
    } else {
      censored.ic50[[ab.tmp]][seqname.index] <- NA
      imputed.ic50[[ab.tmp]][seqname.index] <- NA
    }
    if (length (data.assay.tmp[data.assay.tmp[, 2] != "", 2])) {

      if (">" %in% substr (data.assay.tmp[, 2], 1, 1)) {
        censored.ic80[[ab.tmp]][seqname.index] <- 1
      }
      imputed.ic80[[ab.tmp]][seqname.index] <- merge.readouts (data.assay.tmp[, 2])
    } else {
      censored.ic80[[ab.tmp]][seqname.index] <- NA
      imputed.ic80[[ab.tmp]][seqname.index] <- NA
    }
  }
}

# combine results into new dataframe
readouts <- data.frame (seq.id.catnap=seqname.selected.db)
for (ab.tmp in antibodies) {
  readouts.tmp <- data.frame (censored.ic50[[ab.tmp]], censored.ic80[[ab.tmp]],
                              imputed.ic50[[ab.tmp]], imputed.ic80[[ab.tmp]])
  names (readouts.tmp) <- c (paste0 (ab.tmp, ".ic50.censored"),
                             paste0 (ab.tmp, ".ic80.censored"),
                             paste0 (ab.tmp, ".ic50.imputed"),
                             paste0 (ab.tmp, ".ic80.imputed"))
  readouts <- data.frame (readouts, readouts.tmp)
}

# use additive method to determine predicted combinations of IC50/IC80 
# (i.e., the "Quantitative 1" endpoint)
readouts$pc.ic50 <- apply (readouts[, grep ("ic50.imputed", names (readouts), fixed = TRUE)], 1, wagh.additive.method)
readouts$log10.pc.ic50 <- log10 (readouts$pc.ic50)
readouts$pc.ic80 <- apply (readouts[, grep ("ic80.imputed", names (readouts), fixed = TRUE)], 1, wagh.additive.method)
readouts$log10.pc.ic80 <- log10 (readouts$pc.ic80)

# calculate IIP (a.k.a. the "Quantitive 2" endpoint)
iip.c <- 10
iip.m <- log10 (4) / (readouts$log10.pc.ic80 - readouts$log10.pc.ic50)
iip.f.c <- (iip.c ^ iip.m) / ((readouts$pc.ic50 ^ iip.m) + (iip.c ^ iip.m))
iip.f.c[iip.f.c >= 1] <- 1 - .Machine$double.neg.eps
readouts$iip <- -log10 (1 - iip.f.c)

# derive the "Dichotomous 1" endpoint (i.e., is the PC IC50 higher than the
# sensitivity cutoff?)
sensitivity.threshold <- 1
readouts$dichotomous.1 <- as.numeric (readouts$pc.ic50 >= sensitivity.threshold)

# derive the "Dichotomous 2" endpoint (i.e., is the imputed IC50 greater than 
# the sensitivity threshold for at least two antibodies?) (the use of "two" 
# Abs is used whether or not we have a two-Ab cocktail, or three (or more)-Ab
# cocktail, etc., as per earlier discussion)
min.resistant.abs <- 2
readouts$dichotomous.2 <- as.numeric (rowSums (readouts[, grep ("ic50.imputed", names (readouts), fixed = TRUE)] >= sensitivity.threshold) >= min.resistant.abs)


# ---------------------------------------------------------------------------- #
# STEP 2:  take care of business
# ---------------------------------------------------------------------------- #


# let's start off by making our HXB2 map
hxb2.seq <- toupper (paste (data.seq$B.FR.1983.HXB2.K03455, collapse="")) # does not seem like robust mapping
hxb2.map <- mk.hxb2.map (hxb2.seq)

# generate our AA residue information and remove invariant columns
env.aa.indicators <- as.data.frame (aa.char.2.aa.binary (seqs.selected, seqname.selected.full, hxb2.map))
env.aa.indicators.reduced <- env.aa.indicators[, !(colMeans (env.aa.indicators) %in% 0:1)]

# assemble our data set
data.final <- data.frame (seq.id.lanl = seqname.selected.full,
                          seq.id.catnap = seqname.selected.db,
                          readouts[, -1],
                          as.data.frame (country.2.geo.reg (country.selected)),
                          as.data.frame (bin.subtype (subtype.selected)),
                          env.aa.indicators.reduced)

# name our outfile
filename <- paste0 ("multiab_catnap_", paste (antibodies, collapse = "_"), "_", format (Sys.time (), "%d%b%Y"), ".csv")
setwd (path.data.analysis)
write.csv (data.final, file = filename, row.names = FALSE)


# ---------------------------------------------------------------------------- #
#                                    - 30 -
# ---------------------------------------------------------------------------- #