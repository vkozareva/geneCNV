#!/bin/bash
# Train the model on a set of normal bam files.
(cd ..; python setup.py install)
cnv create-matrix training_bam.fofn --wanted_gene DMD --outfile training_coverage_matrix_wanted.csv --targetfile training_targets_wanted.pickle
cnv create-matrix training_bam.fofn --targets_bed_file training_targets.bed --outfile training_coverage_matrix.csv --targetfile training_targets.pickle
diff training_coverage_matrix_wanted.csv training_coverage_matrix.csv
diff training_targets_wanted.pickle training_targets.pickle
cnv train-model training_targets.pickle training_coverage_matrix.csv training_parameters.pickle
