"""The main Metropolis Hastings Model """

import logging
import numpy as np
import scipy
from scipy.signal import savgol_filter

from IntensitiesDistribution import IntensitiesDistribution
from TargetJointDistribution import TargetJointDistribution
from CopyNumberDistribution import CopyNumberDistribution

class PloidyModel(object):
    """This is the full statistical model and class that runs the Metropolis Hastings sampling scheme. It is responsible for taking a
    parameter set, subject data, and copy number support and running MCMC to determine the
    posterior probability of different ploidy states.
    cnv_support -- array-like containing the different possible ploidy states (ints)
    hln_parameters -- instance of HLN_Parameters containing mu (array), covariance (matrix), and targets(list)
    """

    def __init__(self, cnv_support, hln_parameters, data=None, ploidy=None, intensities=None, first_baseline_i=None,
                 exclude_covar=False):
        """Initialize the data model with its input arguments.
        Load the parameters and initialize starting states as necessary."""

        self.mu = hln_parameters.mu
        self.covariance = hln_parameters.covariance
        self.targets = hln_parameters.targets
        self.n_targets = len(self.targets)
        self.first_baseline_i = self.n_targets if first_baseline_i is None else first_baseline_i
        self.data = data

        self.cnv_support = cnv_support

        # initialize values
        self.initStates(ploidy, intensities)
        self.likelihoods = None

        # initialize joint distribution with data and parameters
        self.joint_target = TargetJointDistribution(self.mu, self.covariance, self.cnv_support, self.data,
                                                    exclude_covar=exclude_covar)

    def initStates(self, ploidy=None, intensities=None):
        """Reset the ploidy and intensity states"""
        self.intensities = IntensitiesDistribution(self.mu, self.covariance).sample() if intensities is None else intensities
        self.ploidy = CopyNumberDistribution(self.n_targets,
                                             support=self.cnv_support).sample_prior(self.first_baseline_i) if ploidy is None else ploidy

    def RunMCMC(self, n_iterations=10000, prior_copy_data=None, prior_mcmc_intens=None, prior_likelihoods=None):
        """Metropolis Hastings sampling of the posterior likelihood"""
        self.mcmc_copy_data = np.zeros((self.n_targets, n_iterations))
        self.mcmc_intens = np.zeros((n_iterations, self.n_targets))
        self.likelihoods = np.zeros(n_iterations)
        self.acceptance = np.zeros((n_iterations, self.n_targets))

        # Targets with labels beginning with 'Baseline' only have their intensities sampled.
        for i in xrange(n_iterations):
            for target_i in xrange(self.n_targets):
                self.ploidy[target_i], self.intensities[target_i], self.acceptance[i, target_i] = self.joint_target.sample(
                    self.ploidy, self.intensities, target_i, target_i >= self.first_baseline_i)
                self.mcmc_copy_data[target_i, i] = self.ploidy[target_i]
                self.mcmc_intens[i, target_i] = self.intensities[target_i]

            self.likelihoods[i] = self.joint_target.log_joint_likelihood(self.intensities, self.ploidy)

            # Log some convergence info at decile intervals.
            if (i + 1) % (n_iterations / 10) == 0:
                logging.info('After {} iterations:\ncnv: {}\nlikelihood: {}\n'.format(
                    i + 1, self.ploidy, self.likelihoods[i]))

        # Log acceptance ratio at end
        logging.info('Acceptance ratio: {}'.format(np.mean(self.acceptance)))

        # Combine with any previously computed sampling data
        # all or none should be passed in
        if prior_copy_data is not None:
            self.mcmc_copy_data = np.concatenate((prior_copy_data, np.copy(self.mcmc_copy_data)), axis=1)
            self.mcmc_intens = np.concatenate((prior_mcmc_intens, np.copy(self.mcmc_intens)), axis=0)
            self.likelihoods = np.concatenate((prior_likelihoods, np.copy(self.likelihoods)))

            logging.info('Using previously passed iteration data, updating to {} total iterations'.format(len(self.likelihoods)))

    def ReportMCMCData(self, burn_in=1000, autocor_slice=100):
        """Report on the posterior distribution obtained by the sampling procedure,
           incorporating burn-in and autocorrelation corrections. """

        self.copy_posteriors = np.zeros((self.n_targets, len(self.cnv_support)))

        for target_i in xrange(self.n_targets):
            # Exclude samples before burn in and then take only every 100th sample to reduce autocorrelation.
            copy_slice = self.mcmc_copy_data[target_i][burn_in:][::autocor_slice]
            for cni, copy_num in enumerate(self.cnv_support):
                self.copy_posteriors[target_i][cni] = np.sum(copy_slice == copy_num)
        self.copy_posteriors /= float(len(copy_slice))

        # find some way to return the intensity distributions?
        return self.copy_posteriors

    def LikelihoodComparison(self, norm_copy_num):
        """Returns the difference in intensity-optimized log-likelihood between the ploidy state comprised of
        individual target modes and the expected normal ploidy state (determined by norm_copy_num). Can only
        be called after ReportMCMCData().

        Larger values indicate greater confidence in the combined mode ploidy state. Negative values large in
        magnitude indicate multimodality in target distribution and that sampler has gotten stuck in
        non-optimal mode.

        norm_copy_num -- The normal ploidy number for all targets in a non-carrier individual
        """
        # get mode ploidy state for each target and generate normal ploidy state
        MAP_ploidy = np.take(self.cnv_support, np.argmax(self.copy_posteriors, axis=1))
        normal_ploidy = norm_copy_num * np.ones(self.n_targets)

        args_map = (MAP_ploidy, True)
        args_norm = (normal_ploidy, True)

        # optimize joint log likelihood for intensities given copy numbers and data
        optarg_map = scipy.optimize.minimize(self.joint_target.log_joint_likelihood, np.concatenate((self.mu.flatten(), [0.])),
                                             args=args_map, tol=1e-6)
        optarg_norm = scipy.optimize.minimize(self.joint_target.log_joint_likelihood, np.concatenate((self.mu.flatten(), [0.])),
                                              args=args_norm, tol=1e-6)

        log_like_diff = (self.joint_target.log_joint_likelihood(optarg_map.x, MAP_ploidy) -
                         self.joint_target.log_joint_likelihood(optarg_norm.x, normal_ploidy))
        return log_like_diff

    def DetectModeJump(self, window_length=201, polyorder=5, initial_offset=500):
        """Detect jumps in log likelihood indicative of switching from one metastable mode to another during iterations.
        Uses Savitzky-Golay filter to first smooth noisy log-likelihood data, then finds highest peak (corresponding to sharpest
        jump in log-likelihood) in numerically differentiated data.
        Method should only be called after sampling has run and if metastability error suspected.

        Returns index of peak and height of peak (in differentiated data). If metastability error not detectable
        will returned peak will likely be due to noise and will have low height.

        window_length -- length of filter window as in scipy.signal.savgol_filter (must be odd, positive int)
        polyorder -- polynomial order as in scipy.signal.savgol_filter (must be smaller than window_length)
        initial_offset -- initial iteration data to remove from peak analysis (similar to burn-in period when likelihood
                          increases significantly)
        """
        smooth_likelihoods = savgol_filter(self.likelihoods, window_length, polyorder)
        grad_smooth = np.gradient(smooth_likelihoods)

        # return index and height of largest peak in gradient
        peak_ind = np.argmax(grad_smooth[initial_offset:]) + initial_offset
        peak_height = np.amax(grad_smooth[initial_offset:])

        return peak_ind, peak_height


