#! /bin/bash

Rscript -e "getwd(); list.files('${TRAVIS_BUILD_DIR}')"
Rscript -e "setwd('${TRAVIS_BUILD_DIR}/docs/'); bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
Rscript -e "setwd('${TRAVIS_BUILD_DIR}/docs/'); bookdown::render_book('index.Rmd', 'bookdown::pdfbook')"