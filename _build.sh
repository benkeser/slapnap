!# /bin/bash

Rscript -e "setwd('${TRAVIS_BUILD_DIR}/docs/'); bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
Rscript -e "setwd('${TRAVIS_BUILD_DIR}/docs/'); bookdown::render_book('index.Rmd', 'bookdown::pdfbook')"