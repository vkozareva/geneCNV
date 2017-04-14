#!/bin/bash
# Train the model on a set of normal bam files.
(cd ..; python setup.py install)
cnv create-matrix training_bam.fofn --wanted_gene DMD --outfile training_coverage_matrix_wanted.csv --targetfile training_targets_wanted.pickle
cnv create-matrix training_bam.fofn --targets_bed_file DMD_targets.bed --outfile DMD_coverage_matrix.csv --targetfile DMD_targets.pickle
cnv create-matrix training_bam.fofn --targets_bed_file DMD_with_baseline_targets.bed --outfile DMD_with_baseline_coverage_matrix.csv --targetfile DMD_with_baseline_targets.pickle
diff training_coverage_matrix_wanted.csv DMD_coverage_matrix.csv
diff training_targets_wanted.pickle DMD_targets.pickle
cnv train-model training_targets_wanted.pickle training_coverage_matrix_wanted.csv training_wanted_parameters.pickle
cnv train-model DMD_targets.pickle DMD_coverage_matrix.csv DMD_parameters.pickle
cnv train-model DMD_with_baseline_targets.pickle DMD_with_baseline_coverage_matrix.csv DMD_with_baseline_parameters.pickle
