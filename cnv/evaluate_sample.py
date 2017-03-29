import sys

import logging
import cPickle
import numpy as np
import pandas as pd
from mando import main
from mando import command

from coverage_matrix import CoverageMatrix
from Gibbs.PloidyModel import PloidyModel
from LogisticNormal import hln_EM
from hln_parameters import HLN_Parameters
from Gibbs.VisualizeMCMC import VisualizeMCMC


#pylint: disable=unused-argument

@command('evaluate-sample')
def evaluate_sample(subjectBamfilePath, parametersFile, outputFile, n_iterations=10000, exclude_covar=False, norm_cutoff=0.3):
    """Test for copy number variation in a given sample

    :param subjectBamfilePath: The ID of the subject
    :param parametersFile: The files with the model parameters
    :param outputFile: Output file name
    :param n_iterations: The number of MCMC iterations desired
    :param exclude_covar: If True, exclude covariance estimates in calculations of conditional and joint probabilities
    :param norm_cutoff: The cutoff for posterior probability of the normal target copy number, below which targets are flagged

    """
    logging.info("Running evaluate samples")

    # Read the parameters file.
    hln_parameters = cPickle.load(open(parametersFile, 'rb'))
    targets = hln_parameters.targets

    # Parse subject bamfile
    subject_df = CoverageMatrix().create_coverage_matrix([subjectBamfilePath], targets)
    subject_id = subject_df['subject'][0]
    # get 'normal' copy number based on whether subject is male or female
    norm_copy_num = 1 if subject_id[0] == 'M' else 2
    target_columns = [target['label'] for target in targets]
    subject_data = subject_df[target_columns].values.astype('float')

    # add option to expand this support later?
    cnv_support = [1, 2, 3]

    ploidy_model = PloidyModel(cnv_support, hln_parameters, data=subject_data, exclude_covar=exclude_covar)
    ploidy_model.RunMCMC(n_iterations)
    copy_posteriors = ploidy_model.ReportMCMCData()

    mcmc_df = pd.DataFrame(copy_posteriors, columns=['copy_{}'.format(cnv) for cnv in cnv_support])
    mcmc_df['Target'] = target_columns
    # do we want to print this dataframe to a csv or output file? -- otherwise I can just remove it

    visualize_instance = VisualizeMCMC(cnv_support, target_columns, copy_posteriors)
    visualize_instance.visualize_copy_numbers('Copy Number Posteriors for Subject {}'.format(subject_id), outputFile)

    for target_i in xrange(len(copy_posteriors)):
        norm_i = cnv_support.index(norm_copy_num)
        if copy_posteriors[target_i][norm_i] < norm_cutoff:
            high_copy = np.argmax(copy_posteriors[target_i])
            high_posterior = copy_posteriors[target_i][high_copy]
            logging.info('{} has a posterior probability of {} of having {} copies'.format(target_columns[target_i], high_posterior, cnv_support[high_copy]))

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
    mu, covariance = hln_EM(np.array(coverage_df[targetCols].values).astype(float), max_iterations=150, tol=1e-11)

    # Pickle the intervals and hyperparameters into the outputFile.
    logging.info('Writing intervals plus hyperparameters to file {}.'.format(outputFile))
    hln_parameters = HLN_Parameters(targets, mu, covariance)
    with open(outputFile, 'w') as f:
        cPickle.dump(hln_parameters, f)

if __name__ ==  '__main__':
    main()
