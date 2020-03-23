#!/bin/sh

# post-build test script

# check options checking
Rscript tests/test-check_opts.R

# check utility functions
Rscript tests/test-utils.R

# check analysis dataset
Rscript tests/test-analysis_dataset.R

# check *small* SL run
Rscript tests/test-run_sl.R

# check *small* vimp run
Rscript tests/test-vimp.R

# ------------------------------------------
# Single nAb, single outcome, single learner
# ------------------------------------------
# no cv, no vimp, ic50

# no cv, no vimp, sens1

# no cv, vimp, ic50

# cv, vimp, ic50

# ------------------------------------------
# Multiple nAbs, single outcome, single learner
# ------------------------------------------
# no cv, no vimp, ic50

# no cv, no vimp, iip

# no cv, vimp

# cv, vimp

# ------------------------------------------
# Multiple nAbs, single outcome, multiple learners
# ------------------------------------------
# no cv, no vimp, ic50

# no cv, no vimp, iip

# no cv, vimp

# cv, vimp
