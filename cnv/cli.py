import sys

import logging
import cPickle
import numpy as np
import pandas as pd
from mando import main
from mando import command

import utilities as cnv_util
from coverage_matrix import CoverageMatrix
from MCMC.PloidyModel import PloidyModel
from LogisticNormal import hln_EM
from hln_parameters import HLN_Parameters
from MCMC.VisualizeMCMC import VisualizeMCMC


@command('create-matrix')
def create_matrix(bamfiles_fofn, outfile=None, target_argfile=None, wanted_gene=None, targets_bed_file=None, unwanted_filters=None, min_dist=629):
    """ Create coverage_matrix from given bamfiles_fofn.

    :param bamfiles_fofn: File containing the paths to all bedfiles to be included in the coverage_matrix
    :param outfile: The path to a csv output file to create from the coverage_matrix. If not provided, no output file will be created.
    :param target_argfile: Path to an output file to contain pickled dict holding target intervals, unwanted_filters, and min_dist.
    :param wanted_gene: Gene from which to gather targets
    :param targets_bed_file: Alternative source of targets, and that may include baseline intervals
    :param unwanted_filters: Comma separated list of filters on reads that should be skipped, keyed by the name of the filter
    :param min_dist: Any two intervals that are closer than this distance will be merged together,
        and any read pairs with insert lengths greater than this distance will be skipped. The default value of 629
        was derived to be one less than the separation between intervals for Exon 69 and Exon 70 of DMD.

    Valid filter names: unmapped, MAPQ_below_60, PCR_duplicate, mate_is_unmapped, not_proper_pair, tandem_pair,
                        negative_insert_length, insert_length_greater_than_merge_distance, pair_end_less_than_reference_end

    """

    if bamfiles_fofn.endswith('.bam'):
        bamfiles_fofn = bamfiles_fofn.split(',')
    if wanted_gene and targets_bed_file:
        logging.error("Both --wanted_gene and --targets_bed_file were specified (only one is allowed).")
        sys.exit(1)
    elif wanted_gene:
        targets = cnv_util.combine_panel_intervals(wanted_gene=wanted_gene, min_dist=min_dist)
    elif targets_bed_file:
        targets = []
        with open(targets_bed_file) as f:
            for line in f:
                line = line.rstrip('\r\n')
                chrom, start, end, name = line.split('\t')
                targets.append({'chrom': chrom, 'start': int(start), 'end': int(end), 'label': name})
    else:
        logging.error("One of --wanted_gene or --targets_bed_file must be specified.")
        sys.exit(1)

    if unwanted_filters is not None:
        unwanted_filters = unwanted_filters.split(',')

    if target_argfile:
        targets_params = {'targets': targets,
                          'unwanted_filters': unwanted_filters,
                          'min_dist': min_dist}
        with open(target_argfile, 'w') as f:
            cPickle.dump(targets_params, f, protocol=cPickle.HIGHEST_PROTOCOL)

    matrix_instance = CoverageMatrix(unwanted_filters=unwanted_filters, min_interval_separation=min_dist)
    coverage_matrix_df = matrix_instance.create_coverage_matrix(bamfiles_fofn, targets)
    if outfile:
        coverage_matrix_df.to_csv(outfile)
        print 'Finished creating {}'.format(outfile)


@command('evaluate-sample')
def evaluate_sample(subjectBamfilePath, parametersFile, outputFile, n_iterations=10000, exclude_covar=False, norm_cutoff=0.3):
    """Test for copy number variation in a given sample

    :param subjectBamfilePath: Path to subject bamfile (.bam.bai must be in same directory)
    :param parametersFile: Pickled file containing a dict with CoverageMatrix arguments and
                           instance of HLN_Parameters (mu, covariance, targets)
    :param outputFile: Output file name without extension -- generates two output files (one .txt file of posteriors
                       and one .pdf displaying stacked bar chart)
    :param n_iterations: The number of MCMC iterations desired
    :param exclude_covar: If True, exclude covariance estimates in calculations of conditional and joint probabilities
    :param norm_cutoff: The cutoff for posterior probability of the normal target copy number, below which targets are flagged

    """
    logging.info("Running evaluate samples")

    # Read the parameters file.
    targets_params = cPickle.load(open(parametersFile, 'rb'))
    targets = targets_params['parameters'].targets

    # Parse subject bamfile
    subject_df = CoverageMatrix(unwanted_filters=targets_params['unwanted_filters'],
                                min_interval_separation=targets_params['min_dist']).create_coverage_matrix([subjectBamfilePath], targets)
    subject_id = subject_df['subject'][0]
    # get 'normal' copy number based on whether subject is male or female
    norm_copy_num = 1. if subject_id[0] == 'M' else 2.
    target_columns = [target['label'] for target in targets]
    subject_data = subject_df[target_columns].values.astype('float')

    # add option to expand this support later?
    # note that having 0 in support causes problems in the joint probability calculation
    cnv_support = [1e-10, 1, 2] if subject_id[0] == 'M' else [1, 2, 3]
    # until normalization against other genes, initializing with most probable normal state
    # These states are fixed for the baseline targets.
    initial_ploidy = norm_copy_num * np.ones(len(targets))

    ploidy_model = PloidyModel(cnv_support, targets_params['parameters'], data=subject_data,
                               ploidy=initial_ploidy, exclude_covar=exclude_covar)
    ploidy_model.RunMCMC(n_iterations)
    copy_posteriors = ploidy_model.ReportMCMCData()

    mcmc_df = pd.DataFrame(copy_posteriors, columns=['Copy_{}'.format(cnv) for cnv in cnv_support])
    mcmc_df['Target'] = target_columns
    mcmc_df.to_csv('{}.txt'.format(outputFile), sep='\t')

    visualize_instance = VisualizeMCMC(cnv_support, target_columns, copy_posteriors)
    visualize_instance.visualize_copy_numbers('Copy Number Posteriors for Subject {}'.format(subject_id), '{}.pdf'.format(outputFile))

    for target_i in xrange(len(copy_posteriors)):
        norm_i = cnv_support.index(norm_copy_num)
        if copy_posteriors[target_i][norm_i] < norm_cutoff:
            high_copy = np.argmax(copy_posteriors[target_i])
            high_posterior = copy_posteriors[target_i][high_copy]
            logging.info('{} has a posterior probability of {} of having {} copies'.format(target_columns[target_i],
                                                                                           high_posterior, cnv_support[high_copy]))


@command('train-model')
def train_model(targetsFile, coverageMatrixFile, outputFile, fit_diag_only=False):
    """Train a model that detects copy number variation.

    :param targetsFile: Pickled file containing target intervals and CoverageMatrix arguments
    :param coverageMatrixFile: CSV file containing coverage data for all samples of interest
    :param outputFile: Output file name, returns CoverageMatrix arguments, and HLN_Parameters object in pickled dict
    :param fit_diag_only: if True, will return diagonal matrix after fitting only variances (all off-diag 0)
    """
    logging.info("Running sample training.")
    if fit_diag_only:
        logging.info('Fitting only diagonal variance terms.')

    # Read the targets file.
    with open(targetsFile) as f:
        targets_params = cPickle.load(f)
        targets = targets_params['targets']

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
                logging.warning('Subject {} has no coverage for target {}.'.format(subject['id'], targetCol))
    if errors > 0:
        sys.exit(1)

    # Compute the logistic normal hyperparameters.
    # Omit the non-target columns.
    mu, covariance = hln_EM(np.array(coverage_df[targetCols].values).astype(float), max_iterations=150, tol=1e-11,
                            fit_diag_only=fit_diag_only)

    # Pickle the intervals, hyperparameters and CoverageMatrix arguments into the outputFile.
    logging.info('Writing intervals plus hyperparameters to file {}.'.format(outputFile))
    targets_params['parameters'] = HLN_Parameters(targets, mu, covariance)
    del targets_params['targets']
    with open(outputFile, 'w') as f:
        cPickle.dump(targets_params, f, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    main()
