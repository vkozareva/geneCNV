import sys

import logging
import cPickle
import numpy as np
import pandas as pd
from mando import main
from mando import command

import utilities as cnv_util
from coverage_matrix import CoverageMatrix
from LogisticNormal import hln_EM
from hln_parameters import HLN_Parameters
from MCMC.VisualizeMCMC import VisualizeMCMC
from MCMC.ConvergenceAnalysis import ConvergenceAnalysis

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
        targets_params = {'full_targets': targets,
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
def evaluate_sample(subjectBamfilePath, parametersFile, outputFile, n_iterations=10000, burn_in_prop=0.3, autocor_slice=50,
                    exclude_covar=False, no_gelman_rubin=False, num_chains=4, max_iterations=25000, norm_cutoff=0.5):
    """Test for copy number variation in a given sample

    :param subjectBamfilePath: Path to subject bamfile (.bam.bai must be in same directory)
    :param parametersFile: Pickled file containing a dict with CoverageMatrix arguments and
                           instance of HLN_Parameters (mu, covariance, targets)
    :param outputFile: Output file name without extension -- generates two output files (one .txt file of posteriors
                       and one .pdf displaying stacked bar chart)
    :param n_iterations: The number of MCMC iterations desired (should be divisible by 100)
    :param burn_in_prop: The proportion of MCMC iterations to exclude as part of burn-in period (should be divisible by 0.05)
    :param autocor_slice: The autocorrelation slice coefficient to use when reporting posterior probabilities
                          ie. only every 50th iteration will be kept
    :param exclude_covar: If True, exclude covariance estimates in calculations of conditional and joint probabilities
    :param no_gelman_rubin: Will not do Gelman-Rubin convergence analysis before metastability analysis
    :param num_chains: Number of independent chains to use during G-R analysis
    :param max_iterations: Maximum number of iterations to use during convergence analysis (both G-R and metastability)
    :param norm_cutoff: The cutoff for posterior probability of the normal target copy number, below which targets are flagged

    """
    logging.info("Running evaluate samples")

    # Read the parameters file.
    targets_params = cPickle.load(open(parametersFile, 'rb'))
    full_targets = targets_params['full_targets']
    targets_to_test = targets_params['parameters'].targets

    # Parse subject bamfile
    subject_df = CoverageMatrix(unwanted_filters=targets_params['unwanted_filters'],
                                min_interval_separation=targets_params['min_dist']).create_coverage_matrix([subjectBamfilePath], full_targets)
    subject_id = subject_df['subject'][0]

    # get 'normal' copy number based on whether subject is male or female
    norm_copy_num = 1. if subject_id[0] == 'M' else 2.
    target_columns = [target['label'] for target in targets_to_test]
    subject_data = subject_df[target_columns].values.astype('float')

    # Note that having 0 in support causes problems in the joint probability calculation if off-target reads exist
    cnv_support = np.array([1e-10, 1, 2, 3]).astype(float)

    # Check if baseline targets exist
    first_baseline_i = len(targets_to_test) # In case there are no baseline targets.
    for i in xrange(len(targets_to_test)):
        if targets_to_test[i]['label'].startswith('Baseline'):
            first_baseline_i = i
            logging.info('Using {} {} baseline targets to normalize ploidy number'.format(('sum of' if 'Sum' in targets_to_test[i]
                                                                                           ['label'] else 'individual'),
                                                                                           len(full_targets) - first_baseline_i))
            break

    # ploidy model (and sampling) actually run within convergence analysis instance
    convergence_analysis = ConvergenceAnalysis(cnv_support, targets_params['parameters'], subject_data, first_baseline_i,
                                               exclude_covar, n_iterations, burn_in_prop)
    if not no_gelman_rubin:
        convergence_analysis.gelman_rubin_analysis(num_chains, len(targets_to_test), max_iterations=max_iterations)

    # Check whether result is far from optimal mode (assuming normal ploidy) and repeat to avoid metastability error
    # note that this will only catch metastabality errors that lead to false positives, not false negatives
    copy_posteriors, loglike_diff = convergence_analysis.metastability_error_analysis(norm_copy_num, autocor_slice,
                                                                                      max_iterations=max_iterations)

    logging.info('Difference in optimized mode and expected ploidy likelihoods is {}'.format(loglike_diff))

    # Create dataframe for reporting copy number posteriors
    mcmc_df = pd.DataFrame(copy_posteriors, columns=['Copy_{}'.format(cnv) for cnv in cnv_support])
    mcmc_df['target'] = target_columns
    mcmc_df['chrom'] = [target['chrom'] for target in targets_to_test]
    mcmc_df['start'] = [target['start'] for target in targets_to_test]
    mcmc_df['end'] = [target['end'] for target in targets_to_test]
    mcmc_df.to_csv('{}.txt'.format(outputFile), sep='\t')

    # Create stacked bar plot and write to pdf
    visualize_instance = VisualizeMCMC(cnv_support, target_columns, copy_posteriors[:first_baseline_i])
    visualize_instance.visualize_copy_numbers('Copy Number Posteriors for Subject {}'.format(subject_id), '{}.pdf'.format(outputFile))

    # Log targets which seem to have abnormal copy numbers
    norm_index = np.where(cnv_support == norm_copy_num)[0][0]

    for target_i in xrange(len(copy_posteriors[:first_baseline_i])):
        if copy_posteriors[target_i][norm_index] < norm_cutoff:
            high_copy = np.argmax(copy_posteriors[target_i])
            high_posterior = copy_posteriors[target_i][high_copy]
            logging.info(('{} has a posterior probability of {} of having '
                          '{} copies'.format(target_columns[target_i], high_posterior, cnv_support[high_copy])))

    # create summary file of all mutations
    MAP_ploidy = np.take(cnv_support, np.argmax(copy_posteriors[:first_baseline_i], axis=1))
    # get indices of ploidy changes
    ploidy_change_i = np.concatenate(([0], np.where(MAP_ploidy[:-1] != MAP_ploidy[1:])[0] + 1))
    reporting_df = pd.DataFrame(columns=['subject', 'mutation', 'targets', 'chrom', 'start', 'end', 'ploidy',
                                         'max_posterior', 'min_posterior', 'mean_posterior', 'loglikelihood_diff'])

    # if no mutations detected
    if np.all(ploidy_change_i == 0) and MAP_ploidy[0] == norm_copy_num:
        data = [subject_id, 'NORM', '{}-{}'.format(target_columns[0], target_columns[first_baseline_i - 1]),
                targets_to_test[0]['chrom'], targets_to_test[0]['start'], targets_to_test[first_baseline_i - 1]['end'],
                norm_copy_num, np.amax(copy_posteriors[:first_baseline_i, norm_index]),
                np.amin(copy_posteriors[:first_baseline_i, norm_index]),
                np.mean(copy_posteriors[:first_baseline_i, norm_index]), loglike_diff]
        reporting_df.loc[0] = data

    # full data for combined mutations
    for i, t_index in enumerate(ploidy_change_i):
        if MAP_ploidy[t_index] != norm_copy_num:
            end_t_index = ploidy_change_i[i + 1] - 1 if i < len(ploidy_change_i) else first_baseline_i - 1
            target_set = (target_columns[t_index] if (i == len(ploidy_change_i)-1 or ploidy_change_i[i + 1] == t_index + 1) else
                          '{}-{}'.format(target_columns[t_index], target_columns[end_t_index]))

            cnv_index = np.where(cnv_support == MAP_ploidy[t_index])[0][0]
            posterior_set = copy_posteriors[t_index:(end_t_index+1), cnv_index]

            data = [subject_id, ('DEL' if MAP_ploidy[t_index] < norm_copy_num else 'DUP'), target_set,
                    targets_to_test[0]['chrom'], targets_to_test[t_index]['start'], targets_to_test[end_t_index]['end'],
                    round(MAP_ploidy[t_index]), np.amax(posterior_set), np.amin(posterior_set),
                    np.mean(posterior_set), loglike_diff]
            reporting_df.loc[i] = data

    reporting_df[['start', 'end', 'ploidy']] = reporting_df[['start', 'end', 'ploidy']].applymap(int)
    reporting_df.to_csv('{}_summary.txt'.format(outputFile), sep='\t')

@command('train-model')
def train_model(targetsFile, coverageMatrixFile, outputFile, use_baseline_sum=False, max_iterations=150, tol=1e-8,
                fit_diag_only=False):
    """Train a model that detects copy number variation.

    :param targetsFile: Pickled file containing target intervals and CoverageMatrix arguments
    :param coverageMatrixFile: CSV file containing coverage data for all samples of interest
    :param outputFile: Output file name, returns CoverageMatrix arguments, and HLN_Parameters object in pickled dict
    :param use_baseline_sum: Train on sum of baseline targets, instead of each baseline target individually, will return error if
                             no baseline targets found
    :param max_iterations: Maximum number of iterations to use during EM routine before termination
    :param tol: Tolerance for convergence at which to terminate during EM routine
    :param fit_diag_only: if True, will return diagonal matrix after fitting only variances (all off-diag 0)
    """
    logging.info("Running sample training.")
    if fit_diag_only:
        logging.info('Fitting only diagonal variance terms.')

    # Read the targets file.
    with open(targetsFile) as f:
        targets_params = cPickle.load(f)
        targets = targets_params['full_targets']

    # Read the coverageMatrixFile.
    coverage_df = pd.read_csv(coverageMatrixFile, header=0, index_col=0)
    # convert dates to datetime
    coverage_df.date_modified = pd.to_datetime(coverage_df.date_modified, unit='s')
    coverage_df['date'] = coverage_df.date_modified.dt.date

    # Log the list of subjects actually used, with ratios, and perhaps some other stats.
    logging.info('Subjects trained:\n{}'.format(coverage_df[['id', 'date', 'gender', 'TSID_ratio']]))

    # Get the appropriate target columns
    targetCols = [target['label'] for target in targets]
    if use_baseline_sum:
        # assuming all targets listed before baselines -- get index of first one
        first_baseline_i = targets.index(next(target for target in targets if target['label'].startswith('Baseline')))
        targetCols = targetCols[:first_baseline_i] + ['BaselineSum']
        targets = targets_params['full_targets'][:first_baseline_i]
        sum_target = {
            'chrom': '{}-{}'.format(targets_params['full_targets'][first_baseline_i]['chrom'],
                                    targets_params['full_targets'][-1]['chrom']),
            'start': None,
            'end': None,
            'label': 'BaselineSum'
        }
        targets.append(sum_target)
    # Run some sanity checks.
    errors = 0
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
    mu, covariance = hln_EM(np.array(coverage_df[targetCols].values).astype(float), max_iterations=max_iterations,
                            tol=tol, fit_diag_only=fit_diag_only)

    # Pickle the intervals, hyperparameters and CoverageMatrix arguments into the outputFile.
    logging.info('Trained for {} total targets'.format(len(targets)))
    logging.info('Writing intervals plus hyperparameters to file {}.'.format(outputFile))
    targets_params['parameters'] = HLN_Parameters(targets, mu, covariance)
    with open(outputFile, 'w') as f:
        cPickle.dump(targets_params, f, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    main()
