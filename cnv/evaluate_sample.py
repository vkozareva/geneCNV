import sys

import logging
import cPickle
import numpy as np
import pandas as pd
from mando import main
from mando import command

from LogisticNormal import hln_EM
from hln_parameters import HLN_Parameters


#pylint: disable=unused-argument

@command('evaluate-sample')
def evaluate_sample(subjectID, parametersFile, outputFile):
    """Test for copy number variation in a given sample

    :param subjectID: The ID of the subjectparametersFile: The files with the model parameters
    :param outputFile: Output file name.

    """
    logging.info("Running evaluate samples")

@command('train-model')
def train_model(targetsFile, coverageMatrixFile, outputFile):
    """Train a model that detects copy number variation.

    :param targets: Pickle file containing target intervals (gene is alternative).
    :param coverageMatrixFile: CSV file containing coverage data for all samples of interest
    :param outputFile: Output file name
    """
    logging.info("Running sample training.")

    # Read the targets file.
    with open(targetsFile) as f:
        targets = cPickle.load(f)

    # Read the coverageMatrixFile.
    coverage_df = pd.read_csv(coverageMatrixFile, header=0, index_col=0)
    # convert dates to datetime
    coverage_df.date_modified = pd.to_datetime(coverage_df.date_modified, unit='s')
    coverage_df['date'] = coverage_df.date_modified.dt.date

    # Log the list of subjects actually used, with ratios, and perhaps some other stats.
    logging.info('Subjects trained:\n{}'.format(coverage_df[['id', 'date', 'gender', 'TSID_ratio']]))

    # Run some sanity checks.
    errors = 0
    # Could use a more formal method of obtaining target names.
    targetCols = [target['label'] for target in targets]
    for index, subject in coverage_df.iterrows():
        # Every subject has coverage for every target.
        for targetCol in targetCols:
            if targetCol not in subject:
                logging.error('Subject {} is missing target {}.'.format(subject['id'], targetCol))
                errors += 1
            elif subject[targetCol] == 0:
                logging.error('Subject {} has no coverage for target {}.'.format(subject['id'], targetCol))
                errors += 1
    if errors > 0:
        sys.exit(1)

    # Compute the logistic normal hyperparameters.
    # Omit the non-target columns.
    mu, covariance = hln_EM(np.array(coverage_df[targetCols].values).astype(float), max_iterations=75, tol=1e-11)

    # Pickle the intervals and hyperparameters into the outputFile.
    logging.info('Writing intervals plus hyperparameters to file {}.'.format(outputFile))
    hln_parameters = HLN_Parameters(targets, mu, covariance)
    with open(outputFile, 'w') as f:
        cPickle.dump(hln_parameters, f)

if __name__ ==  '__main__':
    main()
