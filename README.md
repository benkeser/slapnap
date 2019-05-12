
# SLAPNAP

[![Travis-CI Build
Status](https://travis-ci.com/benkeser/slapnap.svg?branch=master)](https://travis-ci.com/benkeser/slapnap)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
<!-- [![DOI](https://zenodo.org/badge/75324341.svg)](https://zenodo.org/badge/latestdoi/75324341) -->

> Super LeArner Predictions using NAb Panels

**Authors:** [David Benkeser](https://www.github.com/benkeser/), Craig
Magaret, Sohail Nizam, Bhavesh Borate, Brian Williamson, Peter Gilbert

-----

## Description

`slapnap` is a Docker image for to develop cross-validation-based
ensemble predictors of neutralization sensitivity/resistance using HIV
sequences from the [CATNAP database](http://www.hiv.lanl.gov/). The
image provides an automated tool for reading the data from the online
database, compiling analytic data sets, developing prediction models,
and summarizing results.

-----

## Usage

This GitHub repository contains the source code needed to build the
slapnap docker image. However, the repository is also set up for
continuous integration via Travis-CI, with built images found on
[DockerHub](https://cloud.docker.com/u/slapnap/repository/docker/slapnap/slapnap).
See the [Docker
website](https://docs.docker.com/docker-for-windows/install/) for
installation instructions.

From a terminal the image can be downloaded from DockerHub via the
command line.

``` bash
docker pull slapnap/slapnap
```

At run time, the user specifies which NAbs are of interest by setting an
environment variable named `Nab` via the `-e` option of `docker run`.
This variable should be set to a semicolon-separated list of Nabs. A
list of possible Nabs included in the CATNAP database can be found
[here](https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/main.comp).
We will provide a reference to a detailed description of the SLAPNAP
workflow at some point in the future. In the end, an HTML report is
produced summarizing the analysis. This report can be accessed on a
local computer by mounting a local drive to `/home/out/` (the directory
in the Docker container where the report is generated) via the `-v`
option to `docker run`. See the [`docker run` help
page](https://docs.docker.com/engine/reference/run/) for more details.

Here is an example for developing a predictor of sensitivity to a
combination of `VRC07-523-LS`, `PGT121`, and `PGDM1400`.

``` bash
docker run -e Nab="VRC07-523-LS;PGT121;PGDM1400" \ 
           -v /path/to/directory:/home/out \ 
           slapnap/slapnap 
```

If the `docker run` command successfully completes, an HTML report will
appear in `/path/to/directory`.

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/benkeser/slapnap/issues).

-----

## License

© 2019- David Benkeser

The contents of this repository are distributed under the MIT license:

    The MIT License (MIT)
    
    Copyright (c) 2019 David C. Benkeser
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
