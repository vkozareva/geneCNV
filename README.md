# geneCNV
[![Build Status](https://circleci.com/gh/GenePeeks/geneCNV.svg?style=shield&circle-token=41203ba7ace9a56592f8070d62b65d0a45fa334c)](https://circleci.com/gh/GenePeeks/geneCNV)
[![Coverage Status](https://coveralls.io/repos/github/GenePeeks/geneCNV/badge.svg?branch=master)](https://coveralls.io/github/GenePeeks/geneCNV?branch=master)

A command-line program for detection of copy number variants using targeted sequencing data.

This repository is divided into several main sections.

1. The `cnv` folder.  This contains the library code that both trains and runs the MCMC sampling algorithm to determine if samples contain duplications or deletions in this gene.
2. The `test_data` folder which contains sample data used in publication validation experiments.

Read the full documentation at [http://genecnv.readthedocs.io](http://genecnv.readthedocs.io).

## Installation

GeneCNV runs on Python 2.7, which you can install for your OS [here](https://wiki.python.org/moin/BeginnersGuide/Download).

### Dependencies
The package also requires the following Python dependencies, which you can install via `pip`.
* [mando](http://pypi.python.org/pypi/mando)
* [matplotlib](http://matplotlib.org)
* [nose](http://pypi.python.org/pypi/nose/1.3.7)
* [NumPy](http://numpy.org)
* [Pandas](http://pandas.pydata.org/)
* [pysam](http://github.com/pysam-developers/pysam)
* [SciPy](http://www.scipy.org/)


To download and install geneCNV (along with any unsatisfied dependencies):
~~~bash
git clone https://github.com/GenePeeks/geneCNV.git
cd geneCNV
pip install -r requirements.txt
python setup.py install
~~~

### Running Tests
Make sure you've installed properly by running unit tests as follows:

~~~bash
./runtests.sh
~~~
The tests may take a few minutes to complete successfully.

## Command Line Interface


~~~bash
genecnv evaluate-samples -h
~~~

