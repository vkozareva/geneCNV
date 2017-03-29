"""The main Gibbs Sampling Model """

import logging
import numpy as np

from IntensitiesDistribution import IntensitiesDistribution
from CopyNumberDistribution import CopyNumberDistribution
from TargetJointDistribution import TargetJointDistribution

class PloidyModel(object):
    """This is the full statistical model and class that runs the Gibbs Sampling. It is responsible for taking a
    parameter set, a BAM file, and a TargetCollection and running a Gibbs sampler to determine the
    posterior probability of different ploidy states."""

    def __init__(self, cnv_support, hln_parameters, data=None, intensities=None, exclude_covar=False):
        """Initialize the data model with its input arguments.
        Load the parameters and calculate the coverage at each interval in the BAM file."""

        self.mu = hln_parameters.mu
        self.covariance = hln_parameters.covariance
        self.targets = hln_parameters.targets
        self.n_targets = len(self.targets)
        self.data = data

        self.cnv_support = cnv_support

        # initialize values
        self.intensities = IntensitiesDistribution(self.mu, self.covariance).sample() if intensities is None else intensities
        # self.ploidy = CopyNumberDistribution(self.n_targets, support=self.cnv_support).sample_prior()
        # until normalization against other genes, initializing with most probable normal state
        self.ploidy = 2.0 * np.ones(self.n_targets)

        # initialize joint distribution with data and parameters
        self.joint_target = TargetJointDistribution(self.mu, self.covariance, self.cnv_support, self.data,
                                                    exclude_covar=exclude_covar)

    def RunMCMC(self, n_iterations=10000):
        """Metropolis Hastings sampling of the posterior likelihood"""
        self.mcmc_copy_data = np.zeros((self.n_targets, n_iterations))
        self.mcmc_intens = np.zeros((n_iterations, self.n_targets))
        self.likelihoods = np.zeros(n_iterations)
        self.acceptance = np.zeros((n_iterations, self.n_targets))

        for i in xrange(n_iterations):
            for target_i in xrange(self.n_targets):
                self.ploidy[target_i], self.intensities[target_i], self.acceptance[i, target_i] = self.joint_target.sample(self.ploidy,
                                                                                                    self.intensities, target_i)
                self.mcmc_copy_data[target_i, i] = self.ploidy[target_i]
                self.mcmc_intens[i, target_i] = self.intensities[target_i]

            self.likelihoods[i] = self.joint_target.log_joint_likelihood(self.ploidy, self.intensities)

            # Log some convergence info at decile intervals.
            if (i + 1) % (n_iterations / 10) == 0:
                logging.info('After {} iterations:\ncnv: {}\nlikelihood: {}\n'.format(
                    i + 1, self.ploidy, self.likelihoods[i]))

        # Log acceptance ratio at end
        logging.info('Acceptance ratio: {}'.format(np.mean(self.acceptance)))

    def ReportMCMCData(self, burn_in=1000, autocor_slice=100):
        """Report on the posterior distribution obtained by the sampling procedure,
           incorporating burn-in and autocorrelation corrections. """

        copy_posteriors = np.zeros((self.n_targets, len(self.cnv_support)))

        for target_i in xrange(self.n_targets):
            # Exclude samples before burn in and then take only every 100th sample to reduce autocorrelation.
            copy_slice = self.mcmc_copy_data[target_i][burn_in:][::autocor_slice]
            copy_posteriors[target_i] = np.bincount(copy_slice.astype(np.int64),
                                                       minlength=len(self.cnv_support) + 1)[1:]
        copy_posteriors /= float(len(copy_slice))

        # find some way to return the intensity distributions?
        return copy_posteriors
