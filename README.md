
# SLAPNAP

> Super LeArner Predictions using NAb Panels

**Authors:** [David Benkeser](https://www.github.com/benkeser/), Craig
Magaret, Sohail Nizam, Bhavesh Borate, Brian Williamson, Peter Gilbert

[![Build
Status](https://travis-ci.com/benkeser/slapnap.svg?token=WgmsWkd2hyf88ZxhK8bp&branch=master)](https://travis-ci.com/benkeser/slapnap)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
<!-- [![DOI](https://zenodo.org/badge/75324341.svg)](https://zenodo.org/badge/latestdoi/75324341) -->

-----

## Description

`slapnap` is a Docker image for developing cross-validation-based
ensemble predictors of neutralization sensitivity/resistance using HIV
sequences from the [CATNAP database](http://www.hiv.lanl.gov/). The
image provides an automated tool for reading the data from the online
database, compiling analytic data sets, developing prediction models,
and summarizing results.

-----

## Usage

This GitHub repository contains the source code needed to build the
slapnap docker image. The repository is also set up for continuous
integration via Travis-CI, with built images found on
[DockerHub](https://cloud.docker.com/u/slapnap/repository/docker/slapnap/slapnap).
See the [Docker
website](https://docs.docker.com/docker-for-windows/install/) for
installation instructions.

From a terminal the image can be downloaded from DockerHub via the
command line.

``` bash
docker pull slapnap/slapnap
```

At run time, the user specifies which nAbs are of interest by setting an
environment variable named `nab` via the `-e` option of `docker run`.
This variable should be set to a semicolon-separated list of nAbs. A
list of possible Nabs included in the CATNAP database can be found
[here](https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/main.comp).
Other options that can be specified using `-e` include:

  - `outcomes`: a semicolon-separated string of outcomes to include in
    the analysis. Possible values are `"ic50"` (included in default),
    `"ic80"`, `"iip"`, `"sens"` (included in default), `"estsens"`,
    `"multsens"`. If only a single `nab` is specified, use `sens` to
    include a dichotomous endpoint. If multiple `nab`s are specified,
    use `estsens` and/or `multsens`. For detailed definitions of
    outcomes, see the “Data details” section of the [full
    documentation](https://benkeser.github.io/slapnap/6-sec-data.html).
  - `sens_thresh`: a numeric value defining the IC\(_{50}\) threshold
    for defining a sensitive versus resistant pseudovirus (default = 1).
    The dichotomous sensitivity/resistant `outcome`s are defined as the
    indicator that (estimated) IC\(_{50}\) is greater than or equal to
    `sens_thresh`.
  - `multsens_nab`: a numeric value used for defining whether a
    pseudovirus is resistant to a multi-nAb cocktail. The dichotomous
    outcome `multsens` is defined as the indicator that a virus has
    (estimated) \(IC_{50}\) greater than `sens_thresh` for at least
    `multsens_nab` nAbs.
  - `learners`: a semicolon-separated string of machine learning
    algorithms to include in the analysis. Possible values include
    `"rf"` (random forest, default), `"xgboost"` (eXtreme gradient
    boosting), and `"lasso"` (elastic net). If more than one algorithm
    is included, then it is assumed that a cross-validated-based
    ensemble (i.e., a super learner) is desired (see the “Method
    details” section of the [full
    documentation](https://benkeser.github.io/slapnap/7-sec-methods.html)).
  - `nfolds`: a numeric string indicating the number of folds to use in
    cross-validation procedures (default = `"2"`).
  - `importance_grp`: a semicolon-separated list of group-level
    biological importance measures to consider (options are none `""`
    \[default\], marginal `"marg"`, conditional `"cond"`, and both). For
    more detail, see the “Variable importance details” section of the
    [full
    documentation](https://benkeser.github.io/slapnap/7-sec-methods.html#variable-importance-details).
  - `importance_ind`: a semicolon-separated list of individual
    variable-level importance measures to consider (options are none
    `""` \[default\], learner-level `"pred"`, and biological marginal
    `"marg"` and conditional `"cond"`, or any combination)
  - `return`: a semicolon-separated list of the output to save in
    addition to the report (options are `"report"` \[default\],
    `"learner"` for the trained algorithm, `"data"` for the analysis
    dataset, `"figures"` for all figures from the report, and `"vimp"`
    for variable importance objects)

For a complete list of options, see the `Dockerfile` and the [full
documentation](https://benkeser.github.io/slapnap/3-sec-runningcontainer.html).

For a detailed description of the SLAPNAP workflow, see the [full
documentation](https://benkeser.github.io/slapnap/).

In the end, an HTML report is produced summarizing the analysis. This
report can be accessed on a local computer by mounting a local drive to
`/home/out/` (the directory in the Docker container where the report is
generated) via the `-v` option to `docker run`. See the [`docker run`
help page](https://docs.docker.com/engine/reference/run/) for more
details.

Here is an example for developing a predictor of sensitivity to a
combination of `VRC07-523-LS`, `PGT121`, and `PGDM1400`.

``` bash
docker run -e nab="VRC07-523-LS;PGT121;PGDM1400" \
           -v /path/to/directory:/home/out \
           slapnap/slapnap
```

If the `docker run` command successfully completes, an HTML report will
appear in `/path/to/directory`. Any errors or messages will also appear
in a nAb combination-specific `.log` file in `/path/to/directory`.

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
