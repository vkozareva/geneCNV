import cPickle
import logging
import sys

import numpy as np
import pandas as pd
from mando import command
from mando import main

from LogisticNormal import hln_EM
from MCMC.ConvergenceAnalysis import ConvergenceAnalysis
from MCMC.VisualizeMCMC import VisualizeMCMC
from cnv import __version__
from cnv.Targets.TargetCollection import DEFAULT_MERGE_DISTANCE
from cnv.Targets.TargetCollection import TargetCollection
from cnv.Targets.Target import Target
from cnv.utilities import SimulateData
from coverage_matrix import CoverageMatrix
from hln_parameters import HLN_Parameters

def configure_logging(verbose=0):
    """ Configure logging and verbosity """
    level = logging.DEBUG if verbose == 2 else logging.INFO if verbose == 1 else logging.WARNING
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=level)

@command('version')
def version():
    """Provide the current version"""
    version_id = 'cnv v{}\n'.format(__version__)
    sys.stdout.write(version_id)
    sys.exit()

@command('create-matrix')
def create_matrix(targetsBedfile, bamfilesFofn, outputFile, targetArgfile=None, unwanted_filters=None,
                  min_dist=DEFAULT_MERGE_DISTANCE, verbose=0):
    """ Create coverage_matrix from given bamfilesFofn.

    :param targetsBedfile: Source of targets, and that may include baseline intervals
    :param bamfilesFofn: File containing the paths to all BAM files to be included in the coverage_matrix
    :param outputFile: The path to a csv output file to create from the coverage_matrix.
    :param targetArgfile: Path to an output file to contain pickled dict holding target intervals, unwanted_filters, and min_dist.
    :param wanted_gene: Gene from which to gather targets
    :param unwanted_filters: Comma separated list of filters on reads that should be skipped, keyed by the name of the filter
    :param min_dist: Any two intervals that are closer than this distance will be merged together,
        and any read pairs with insert lengths greater than this distance will be skipped. The default value of 629
        was derived to be one less than the separation between intervals for Exon 69 and Exon 70 of DMD.
    :param -v, --verbose: 0 - Logging level warning; 1 - Logging level info; 2 - Logging level debug [0]

    Valid filter names: unmapped, MAPQ_below_60, PCR_duplicate, mate_is_unmapped, not_proper_pair, tandem_pair,
                        negative_insert_length, insert_length_greater_than_merge_distance, pair_end_less_than_reference_end

    """
    # set appropriate logging level
    configure_logging(verbose)

    if bamfilesFofn.endswith('.bam'):
        bamfilesFofn = bamfilesFofn.split(',')

    targets = TargetCollection.load_from_txt_file(targetsBedfile, min_merge_dist=min_dist)

    if unwanted_filters is not None:
        unwanted_filters = unwanted_filters.split(',')

    if targetArgfile:
        targets_params = {'full_targets': targets,
                          'unwanted_filters': unwanted_filters
                         }
        with open(targetArgfile, 'w') as f:
            cPickle.dump(targets_params, f, protocol=cPickle.HIGHEST_PROTOCOL)

    matrix_instance = CoverageMatrix(unwanted_filters=unwanted_filters)
    coverage_matrix_df = matrix_instance.create_coverage_matrix(bamfilesFofn, targets)

    coverage_matrix_df.to_csv(outputFile)
    logging.info('Finished creating {}'.format(outputFile))


@command('evaluate-sample')
def evaluate_sample(subjectFilePath, parametersFile, outputPrefix, n_iterations=10000, burn_in_prop=0.3, autocor_slice=50,
                    exclude_covar=False, no_gelman_rubin=False, num_chains=4, use_single_process=False, max_iterations=25000,
                    threshold_loglike_diff=-30, norm_cutoff=0.5, verbose=0):
    """Test for copy number variation in a given sample

    :param subjectFilePath: Path to subject bam (.bam.bai must be in same directory) or coverage count matrix
                            (in csv format) (targets must match those in parametersFile)
    :param parametersFile: Pickled file containing a dict with CoverageMatrix arguments and
                           instance of HLN_Parameters (mu, covariance, targets)
    :param outputPrefix: Output file name without extension -- generates three output files (.txt
                        file of posteriors, _summary.txt, and .pdf with stacked bar chart)
    :param n_iterations: The number of MCMC iterations desired (should be divisible by 100) [10000]
    :param burn_in_prop: The proportion of MCMC iterations to exclude as part of burn-in period
                         (should be divisible by 0.05) [0.3]
    :param autocor_slice: The autocorrelation slice coefficient to use when reporting posterior probabilities
                          ie. only every 50th iteration will be kept [50]
    :param exclude_covar: Exclude covariance estimates in calculations of conditional and joint probabilities
    :param no_gelman_rubin: Will not perform Gelman-Rubin convergence analysis before metastability analysis
    :param num_chains: Number of independent chains to use during G-R analysis, will use separate process for each unless
                       --use_single_process specified [4]
    :param use_single_process: Will not use parallelization during G-R analysis
    :param max_iterations: Maximum number of iterations to use during convergence analysis (both G-R and metastability)
                           [25000]
    :param threshold_loglike_diff: Threshold for calling metastability error in log-likelihood comparison
                                   with normal ploidy state [-30]
    :param norm_cutoff: The cutoff for posterior probability of the normal target copy number, below
                        which targets are flagged [0.5]
    :param -v, --verbose: 0 - Logging level warning; 1 - Logging level info; 2 - Logging level debug [0]

    """
    # set appropriate logging level
    configure_logging(verbose)

    logging.info("Running evaluate samples")

    # Read the parameters file.
    targets_params = cPickle.load(open(parametersFile, 'rb'))
    full_targets = targets_params['full_targets']
    targets_to_test = targets_params['parameters'].targets

    # Parse subject file
    if subjectFilePath.endswith('.csv'):
        subject_df = pd.read_csv(subjectFilePath, index_col=0)
    else:
        subject_df = CoverageMatrix(unwanted_filters=targets_params['unwanted_filters']).create_coverage_matrix([subjectFilePath], full_targets)
    subject_id = subject_df['sample'][0]
    if len(subject_df) > 1:
        logging.warning('Multiple samples in provided CSV. Evaluating only first sample {}.'.format(subject_id))
    target_columns = [target.label for target in targets_to_test]
    # evaluate only first subject if multiple samples in provided CSV
    subject_data = subject_df.iloc[0][target_columns].values.astype('float').flatten()

    # Note that having 0 in support causes problems in the joint probability calculation if off-target reads exist
    cnv_support = np.array([1e-10, 1, 2, 3]).astype(float)

    # Check if baseline targets exist
    first_baseline_i = len(targets_to_test) # In case there are no baseline targets.
    for i in xrange(len(targets_to_test)):
        if 'Baseline' in targets_to_test[i].label:
            first_baseline_i = i
            logging.info('\nUsing {} {} baseline targets to normalize ploidy number'.format(('sum of' if 'Sum' in
                                                                                           targets_to_test[i].label else 'individual'),
                                                                                           len(full_targets) - first_baseline_i))
            break
    # Report non-baseline target coverage
    logging.info('Non-baseline target coverage: {}\n'.format(np.sum(subject_data[:first_baseline_i])))

    # Check and report target-baseline-sum z-score if appropriate
    control_tb_mean = targets_params.get('target_base_mean', None)
    if control_tb_mean:
        control_tb_sd = targets_params.get('target_base_sd', None)
        sample_tb_ratio = np.sum(subject_data[:first_baseline_i]) / subject_data[first_baseline_i]
        sample_zscore = (sample_tb_ratio - control_tb_mean) / control_tb_sd

        logging.info('Z-score of sample target-to-baseline-sum ratio is {}.\n'
                     'Results likely to be less accurate if |z-score| > 1.5.'.format(sample_zscore))
        # set appropriate slice index for correlation check below
        baseline_sum_i = first_baseline_i
    else:
        baseline_sum_i = len(targets_to_test)

    # Check and report test sample and reference set correlation (for non-baseline-sum targets)
    target_xvals = np.exp(np.concatenate((targets_params['parameters'].mu.flatten(), [0]))[:baseline_sum_i])
    control_intensities = target_xvals / np.sum(target_xvals)

    sample_intensities = subject_data[:baseline_sum_i] / np.sum(subject_data[:baseline_sum_i])
    correlation = np.corrcoef(control_intensities, sample_intensities)[0,1]
    logging.info('Correlation between sample intensities and '
                 'average training intensities: {}'.format(round(correlation, 4)))
    if correlation < 0.9:
        logging.warning('Low correlation between test and training samples.\n'
                        'Results likely to be inaccurate if correlation < 0.9.')

    # ploidy model (and sampling) actually run within convergence analysis instance
    convergence_analysis = ConvergenceAnalysis(cnv_support, targets_params['parameters'], subject_data, first_baseline_i,
                                               exclude_covar, n_iterations, burn_in_prop, use_single_process)
    if not no_gelman_rubin:
        convergence_analysis.gelman_rubin_analysis(num_chains, len(targets_to_test), max_iterations=max_iterations)

    # Check whether result is far from optimal mode (assuming normal ploidy) and repeat to avoid metastability error
    # note that this will only catch metastabality errors that lead to false positives, not false negatives
    copy_posteriors, loglike_diff = convergence_analysis.metastability_error_analysis(thresh_loglike_diff=threshold_loglike_diff,
                                                                                      autocor_slice=autocor_slice,
                                                                                      max_iterations=max_iterations)
    norm_copy_num = convergence_analysis.norm_copy_num
    logging.info('Evaluating with normal copy number: {}'.format(norm_copy_num))

    logging.info('Difference in optimized mode and expected ploidy likelihoods is {}'.format(loglike_diff))

    # Create dataframe for reporting copy number posteriors
    mcmc_df = pd.DataFrame(copy_posteriors, columns=['Copy_{}'.format(cnv) for cnv in cnv_support])
    mcmc_df['target'] = target_columns
    mcmc_df['chrom'] = [target.chrom for target in targets_to_test]
    mcmc_df['start'] = [target.start for target in targets_to_test]
    mcmc_df['end'] = [target.end for target in targets_to_test]
    mcmc_df.to_csv('{}.txt'.format(outputPrefix), sep='\t')

    # Create stacked bar plot and write to pdf
    visualize_instance = VisualizeMCMC(cnv_support, target_columns, copy_posteriors[:first_baseline_i])
    visualize_instance.visualize_copy_numbers('Copy Number Posteriors for Subject {}'.format(subject_id), '{}.pdf'.format(outputPrefix))

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
    reporting_df = pd.DataFrame(columns=['subject', 'mutation', 'targets', 'num_targets', 'chrom',
                                         'start', 'end', 'ploidy', 'max_posterior', 'min_posterior',
                                         'mean_posterior', 'loglikelihood_diff', 'name_data'])

    # if no mutations detected
    if np.all(ploidy_change_i == 0) and MAP_ploidy[0] == norm_copy_num:
        posterior_set = copy_posteriors[:first_baseline_i, norm_index]

        name_ind = [i for i,d in enumerate(targets_to_test) if d.name]
        name = '-'.join({targets_to_test[i].name for i in name_ind}) if name_ind else None

        # first_baseline_i also corresponds to length of non-baseline targets in targets_to_test
        data = [subject_id, 'NORM', '{}-{}'.format(target_columns[0], target_columns[first_baseline_i - 1]),
                first_baseline_i, targets_to_test[0].chrom, targets_to_test[0].start,
                targets_to_test[first_baseline_i - 1].end, norm_copy_num,
                np.amax(posterior_set), np.amin(posterior_set), np.mean(posterior_set),
                loglike_diff, name]
        reporting_df.loc[0] = data
    else:
        # full data for combined mutations
        for i, t_index in enumerate(ploidy_change_i):
            if MAP_ploidy[t_index] != norm_copy_num:
                end_t_index = ploidy_change_i[i + 1] - 1 if i < len(ploidy_change_i) - 1 else first_baseline_i - 1
                target_set = (target_columns[t_index] if t_index == end_t_index else
                              '{}-{}'.format(target_columns[t_index], target_columns[end_t_index]))
                num_targets = end_t_index - t_index + 1

                cnv_index = np.where(cnv_support == MAP_ploidy[t_index])[0][0]
                posterior_set = copy_posteriors[t_index:(end_t_index+1), cnv_index]

                # Include target names if existing
                name_ind = [j for j,d in enumerate(targets_to_test[t_index:(end_t_index + 1)]) if d.name]
                name = '-'.join({targets_to_test[j].name for j in name_ind}) if name_ind else None

                data = [subject_id, ('DEL' if MAP_ploidy[t_index] < norm_copy_num else 'DUP'), target_set, num_targets,
                        targets_to_test[0].chrom, targets_to_test[t_index].start, targets_to_test[end_t_index].end,
                        round(MAP_ploidy[t_index]), np.amax(posterior_set), np.amin(posterior_set),
                        np.mean(posterior_set), loglike_diff, name]
                reporting_df.loc[i] = data

    reporting_df[['num_targets', 'start', 'end', 'ploidy']] = reporting_df[['num_targets', 'start',
                                                                            'end', 'ploidy']].applymap(int)

    with open('{}_summary.txt'.format(outputPrefix), 'w') as outfile_main:
        ## outfile_main.write('###### Metadata here\n')
        reporting_df.to_csv(outfile_main, index=False, sep='\t')

@command('train-model')
def train_model(targetsFile, coverageMatrixFile, outputFile, use_baseline_sum=False, max_iterations=150, tol=1e-8,
                fit_diag_only=False, verbose=0):
    """Train a model that detects copy number variation.

    :param targetsFile: Pickled file containing target intervals and CoverageMatrix arguments
    :param coverageMatrixFile: CSV file containing coverage data for all samples of interest
    :param outputFile: Output file name, returns CoverageMatrix arguments, and HLN_Parameters object in pickled dict
    :param use_baseline_sum: Train on sum of baseline targets, instead of each baseline target individually, will return error if
                             no baseline targets found
    :param max_iterations: Maximum number of iterations to use during EM routine before termination [150]
    :param tol: Tolerance for convergence at which to terminate during EM routine [1e-8]
    :param fit_diag_only: Returns diagonal matrix after fitting only variances (all off-diag 0)
    :param -v, --verbose: 0 - Logging level warning; 1 - Logging level info; 2 - Logging level debug [0]
    """
    # set appropriate logging level
    configure_logging(verbose)
    logging.info("Running sample training.")
    if fit_diag_only:
        logging.info('Fitting only diagonal variance terms.')

    # Read the targets file.
    with open(targetsFile) as f:
        targets_params = cPickle.load(f)
        targets = targets_params['full_targets']

    # Read the coverageMatrixFile.
    coverage_df = pd.read_csv(coverageMatrixFile, header=0, index_col=0)

    # Get the appropriate target columns
    targetCols = [target.label for target in targets]
    if use_baseline_sum:
        # assuming all targets listed before baselines -- get index of first one
        first_baseline_i = targets.t_index(next(target for target in targets if 'Baseline' in target.label))
        targetCols = targetCols[:first_baseline_i] + ['BaselineSum']
        # we can slice the TargetCollection
        targets = targets_params['full_targets'][:first_baseline_i]

        chrom_span = '{}-{}'.format(targets_params['full_targets'][first_baseline_i].chrom,
                                    targets_params['full_targets'][-1].chrom)
        sum_target = Target(chrom_span, None, None, 'BaselineSum')
        targets.append(sum_target)

        # Include info about distribution of target:baseline_sum coverage across samples
        target_sums = np.sum(coverage_df[targetCols[:first_baseline_i]].values.astype(float), axis=1)
        target_base_ratios = target_sums.flatten() / coverage_df['BaselineSum'].values.flatten()
        target_base_mean = np.mean(target_base_ratios)
        target_base_sd = np.std(target_base_ratios)
        logging.info('Coefficient of variation of '
                     'target-to-baseline-sum ratio: {}'.format(target_base_sd / target_base_mean))
        targets_params['target_base_mean'] = target_base_mean
        targets_params['target_base_sd'] = target_base_sd

    # Run some sanity checks.
    errors = 0
    for index, subject in coverage_df.iterrows():
        # Every subject has coverage for every target.
        for targetCol in targetCols:
            if targetCol not in subject:
                logging.error('Sample {} is missing target {}.'.format(subject['sample'], targetCol))
                errors += 1
            elif subject[targetCol] == 0:
                logging.warning('Sample {} has no coverage for target {}.'.format(subject['sample'], targetCol))
    if errors > 0:
        sys.exit(1)

    # Compute the logistic normal hyperparameters.
    # Omit the non-target columns.
    mu, covariance = hln_EM(coverage_df[targetCols].values.astype(float), max_iterations=max_iterations,
                            tol=tol, fit_diag_only=fit_diag_only)

    # Pickle the intervals, hyperparameters and CoverageMatrix arguments into the outputFile.
    logging.info('Trained for {} total targets'.format(len(targets)))
    logging.info('Writing intervals plus hyperparameters to file {}.'.format(outputFile))
    targets_params['parameters'] = HLN_Parameters(targets, mu, covariance)

    with open(outputFile, 'w') as f:
        cPickle.dump(targets_params, f, protocol=cPickle.HIGHEST_PROTOCOL)

@command('create-bams')
def create_bams(targetsFile, outputPrefix):
    """Makes simulated data to run the program with, given a target bed file and an output file prefix.

    :param targetsFile: A BED file with targets to simulate coverage for.
    :param outputPrefix: an output file prefix used to name the output bam and fofn file.
    :return:
    """
    SimulateData.make_simulated_data(outputPrefix, targetsFile)


if __name__ == '__main__':
    main()
