# dmd [![CircleCI](https://circleci.com/gh/GenePeeks/dmd.svg?style=svg&circle-token=41203ba7ace9a56592f8070d62b65d0a45fa334c)](https://circleci.com/gh/GenePeeks/dmd)
##The Copy Number Variation Package

A repository dedicated to the analysis of the exons that affect Duchennne and Becker muscular dystrophy. 

This repository is divided into several main sections.

1. The `dmd` folder.  This contains the library code that both trains and runs the Gibbs sampling algorithm to determine if samples contain duplications or deletions in this gene.
2. The `notebooks` folder which contains the analysis scripts and packages we used to develop the algorithm and code base.
3. The `exon_data` folder, which contains sample data used in the notebooks along with graphics developed during the modeling process.

## Installation

~~~bash
git clone https://github.com/GenePeeks/dmd.git
cd dmd
python setup.py install
~~~

## Command Line Interface

~~~bash
dmd evaluate-samples -h
~~~

## Running Tests
~~~bash
./runtests.sh
~~~