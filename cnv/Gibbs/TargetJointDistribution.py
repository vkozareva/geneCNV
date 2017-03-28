import numpy as np
import logging
from IntensitiesDistribution import IntensitiesDistribution
from CopyNumberDistribution import CopyNumberDistribution

class TargetJointDistribution(object):
    """Describes the joint distribution for hierarchical logistic normal model (with multinomial draws).
       Includes methods for calculating unnormalized log likelihood of joint distribution given data
       and sampling both intensity and ploidy for single target conditional on other targets. """

    def __init__(self, targets, mu, covariance, data=None, support=None, exclude_covar=False):
        """Initialize the data to hold past and current samples of ploidy levels for all values in targets."""
        # isinstance(targets, TargetCollection)
        self.data = data
        if support is not None:
            self.support = support
        else:
            self.support = [1, 2, 3]

        self.exclude_covar = exclude_covar

        # compute all inverted components for conditional covariance one time upfront
        if not self.exclude_covar:
            self.cond_inverses = np.zeros((len(mu), len(mu)-1, len(mu)-1))
            for index in xrange(len(mu)):
                cov_22_t = np.concatenate((covariance[:index, :index], covariance[:index, (index+1):]), axis=1)
                cov_22_b = np.concatenate((covariance[(index+1):, :index], covariance[(index+1):, (index+1):]), axis=1)
                cov_22 = np.concatenate((cov_22_t, cov_22_b), axis=0)
                cov_22_inv = np.linalg.inv(cov_22)
                self.cond_inverses[index] = cov_22_inv
            logging.info('Finished computing inverse conditional matrices')

        self.mu = mu
        # keep only diagonal if excluding covariances
        self.covariance = np.diag(np.diagonal(covariance)) if self.exclude_covar else covariance
        self.mu_full = np.concatenate((mu.flatten(), [0]))
        self.inv_covariance_full = np.concatenate((np.concatenate((np.linalg.inv(self.covariance), np.zeros((1,len(self.covariance)))), axis=0),
                                                np.zeros((len(self.covariance)+1,1))), axis=1)

    def sample(self, copies, intensities, target_index):
        """Given a current set of intensities, and the current ploidy state maintained in this class,
        sample a new ploidy state and intensity for a particular target using the prior distribution as the proposal distribution."""
        intensities_proposed = np.copy(intensities)
        copies_proposed = np.copy(copies)

        copy_proposed = CopyNumberDistribution(1, self.support).sample_prior()
        copies_proposed[target_index] = copy_proposed

        # don't update intensity for last target (for identifiability)
        if target_index != (len(self.mu)):
            if self.exclude_covar:
                mu_bar = self.mu[target_index]
                cov_bar = self.covariance[target_index, target_index]
            else:
                mu_bar, cov_bar = self.get_conditional_mvn(self.mu, self.covariance, target_index,
                                                           intensities[:-1], self.cond_inverses[target_index])

            # sample intensity from conditional normal
            intensity_proposed = np.random.normal(mu_bar, np.sqrt(cov_bar))
            intensities_proposed[target_index] = intensity_proposed

            # calculate log likelihood of proposed jump (J(proposed|current state))
            jump_proposed = -(intensity_proposed - mu_bar) ** 2 / (2 * cov_bar)
            jump_previous = -(intensities[target_index] - mu_bar) ** 2 / (2 * cov_bar)
        else:
            intensity_proposed = intensities[target_index]
            jump_proposed = 0
            jump_previous = 0

        # calculate log joint likelihood
        joint_proposed = self.log_joint_likelihood(copies_proposed, intensities_proposed)
        joint_previous = self.log_joint_likelihood(copies, intensities)

        log_test_ratio = joint_proposed + jump_previous - joint_previous - jump_proposed

        # sample for new joint state and keep track of acceptance
        jump_chance = min(1.0, np.exp(log_test_ratio))
        if np.random.rand() < jump_chance:
            return copy_proposed, intensity_proposed, 1.0
        else:
            return copies[target_index], intensities[target_index], 0.0

    def log_joint_likelihood(self, copies, intensities):
        """ Returns unnormalized log likelihood of joint probability given subject data, full set of copy numbers and
            full set of intensities."""
        log_joint = (np.sum(self.data) * -1 * np.log(np.sum(np.multiply(copies, np.exp(intensities)))) +
                     np.sum(np.multiply(self.data, (np.log(copies) + intensities))) +
                     (-0.5 * np.dot(np.dot((intensities - self.mu_full).reshape((1,-1)), self.inv_covariance_full),
                                    (intensities - self.mu_full).reshape((-1,1)))))
        return log_joint

    @staticmethod
    def get_conditional_mvn(mu, cov, index, intensities, cov_22_inv=None):
        """Returns mu and cov for conditional normal distribution for single unknown value
           Pass in previously inverted 2,2 component for faster performance.
           Quadrants are as follows (after translation of desired index to 1,1 position):
           [1,1 [.. ... 1,2 ... ..]]
           [[.. [.. ... ... ... ...]
           [... ... ... ... ... ...]
           [2,1 ... ... 2,2 ... ...]
           [... ... ... ... ... ...]
           [..] ... ... ... ... ..]]

           """
        mu_1 = mu[index]
        mu_2 = np.delete(mu, index)
        a = np.delete(intensities, index)
        # split covariance matrix into quadrants and compute each one separately
        cov_11 = cov[index, index]
        cov_12 = np.concatenate((cov[index][:index], cov[index][(index+1):])).reshape((1, -1))
        cov_21 = np.concatenate((cov[:,index][:index], cov[:, index][(index+1):])).reshape((-1,1))

        if cov_22_inv is None:
            cov_22_t = np.concatenate((cov[:index, :index], cov[:index, (index+1):]), axis=1)
            cov_22_b = np.concatenate((cov[(index+1):, :index], cov[(index+1):, (index+1):]), axis=1)
            cov_22 = np.concatenate((cov_22_t, cov_22_b), axis=0)
            cov_22_inv = np.linalg.inv(cov_22)

        mu_bar = mu_1 + np.dot(np.dot(cov_12, cov_22_inv), (a - mu_2).reshape((-1,1)))
        cov_bar = cov_11 - np.dot(np.dot(cov_12, cov_22_inv), cov_21)

        return mu_bar.flatten()[0], cov_bar.flatten()[0]
