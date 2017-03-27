import sys

import logging
import pickle
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
        targets = pickle.load(f)

    # Read the coverageMatrixFile.
    coverage_df = pd.read_csv(coverageMatrixFile, header=0, index_col=0)
    # convert dates to datetime
    coverage_df.date_modified = pd.to_datetime(coverage_df.date_modified, unit='s')
    coverage_df['date'] = coverage_df.date_modified.dt.date
    # get tsid/tso ratio
    coverage_df['tsid_ratio'] = coverage_df.TSID_only / (coverage_df.TSID_only + coverage_df.TSO_only)

    # Log the list of subjects actually used, with ratios, and perhaps some other stats.
    logging.info('Subjects trained:\n{}'.format(coverage_df[['id', 'date', 'gender', 'tsid_ratio']]))

    # Run some sanity checks.
    errors = 0
    # Could use a more formal method of obtaining target names.
    targetCols = filter(lambda key: key.startswith('Ex'), coverage_df.keys())
    for index, subject in coverage_df.iterrows():
        # Every subject has coverage for every target.
        for targetCol in targetCols:
            if subject[targetCol] == 0:
                logging.error('Subject {} is missing coverage for {}.'.format(subject['id'], targetCol))
                errors += 1
    if errors > 0:
        sys.exit(1)

    # Compute the logistic normal hyperparameters.
    # Omit the non-target columns.
    coverage_df = coverage_df.drop(['tsid_ratio', 'TSID_only', 'TSO_only'], axis=1)
    coverage_grouped_df = coverage_df.groupby('subject').sum()
    mu, covariance = hln_EM(np.array(coverage_grouped_df.values).astype(float), max_iterations=75, tol=1e-11)

    # Pickle the intervals and hyperparameters into the outputFile.
    logging.info('Writing intervals plus hyperparameters to file {}.'.format(outputFile))
    hln_parameters = HLN_Parameters(targets, mu, covariance)
    with open(outputFile, 'w') as f:
        hln_parameters.dump(f)

if __name__ ==  '__main__':
    main()
