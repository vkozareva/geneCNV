#!/bin/bash
# Train the model on a set of normal bam files.
(cd ..; python setup.py install)
coverate_matrix create_matrix training_bam.fofn --outfile training_coverage_matrix.csv --targetfile training_targets.pickle
cnv train-model training_targets.pickle training_coverage_matrix.csv training_parameters.pickle
