# geneCNV
[![Build Status](https://circleci.com/gh/GenePeeks/geneCNV.svg?style=shield&circle-token=41203ba7ace9a56592f8070d62b65d0a45fa334c)](https://circleci.com/gh/GenePeeks/geneCNV)
[![Coverage Status](https://coveralls.io/repos/github/GenePeeks/geneCNV/badge.svg?branch=master)](https://coveralls.io/github/GenePeeks/geneCNV?branch=master)

A command-line program for detection of copy number variants using targeted sequencing data.
GeneCNV is designed for copy number analysis across a subset of genes using parameters
derived from a predefined set of reference (normal) samples.

This repository is divided into several main sections.

1. The `cnv` folder.  This contains the library code that both trains and runs the MCMC sampling algorithm to determine if samples contain duplications or deletions.
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

## Command Line Interface Introduction
GeneCNV involves three main commands: `create-matrix`, `train-model`, and
`evaluate-sample`, corresponding to the following main steps in the
computational pipeline.

### Get coverage counts
To get started, generate coverage counts across relevant targets
and samples using the `create-matrix` command. You must
provide a BED file of relevant targets in this format:
```
X   32867834    32867947    Ex3 DMD
X   33038245    33038327    Ex2 DMD
X   33229388    33229673    Ex1 DMD
```
An example BED file for the DMD gene is provided in `test_data`. Note that the first
four fields (chromosome, start position, end position, label) are required,
while the fifth is optional.

You must also provide a text file of paths to the sample BAM files in this format:
```
/path/to/file1.bam
/path/to/file2.bam
```
An example `create-matrix` command looks like:
~~~bash
genecnv create-matrix test_data/example_dmd_baseline.bed training_samples.fofn \
training_sample_coverage.csv --targetArgfile dmd_baseline_targets.pickle
~~~
Serialized target/argument files can be optionally produced with this command, and
you only need to produce a target/argument file once for a specific set of targets.
An example output CSV for this command is provided in `test_data`. This can be
used to run the subsequent `train-model` command.

### Train the model with normal samples
To estimate the model hyperparameters using all of the samples included in the coverage
count matrix, run the following:
~~~bash
genecnv train-model dmd_baseline_targets.pickle test_data/training_sample_coverage.csv \
dmd_baseline_params.pickle --use_baseline_sum
~~~
Baseline autosomal targets are used to identify absolute copy number when no CNVs are present,
and help provide more accurate results overall. Including baseline targets can also
allow you to identify the sex of a sample when targets on the X chromosome are being
tested. Baseline targets are not analyzed for copy number and are assumed to have
copy number of 2.

If you are using a large number of baseline targets (>20), it's recommended to use
the optional `--use_baseline_sum` argument when calling `train-model`. This
reduces the total number of baseline targets to one during training.

### Evaluate samples for CNVs
Once parameters have been estimated from an appropriate set of training samples,
they can be used to perform copy number analysis for the relevant targets on
a test sample with the `evaluate-sample` command. Here you can pass simply a sample BAM file
or a coverage matrix CSV (generated using the same targets).

To evaluate the first test sample in the file `test_data/test_female_sample_coverage.csv`
use the following command:
~~~bash
genecnv evaluate-sample test_data/test_female_sample_coverage.csv dmd_baseline_params.pickle \
normal_female_results
~~~
This command will produce three output files with the provided prefix, `normal_female_results.txt`,
which provides the posterior probabilities and copy numbers for all relevant targets,
`normal_female_results_summary.txt` which provides a summary of any CNVs detected,
and `normal_female_results.pdf`, which provides a visualization of the copy numbers and
posterior probabilities across targets.

Depending on the number of total targets and MCMC iterations needed for convergence, the
sample evaluation may take up to 10-12 minutes to complete. By default it takes advantage
of multiple cores, but this can be turned off with the option `--use_single_process`.



