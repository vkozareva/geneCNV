import logging
import numpy as np


class TargetJointDistribution(object):
    """Describes the joint distribution for hierarchical logistic normal model (with multinomial draws).
       Includes methods for calculating unnormalized log likelihood of joint distribution given data
       and sampling both intensity and ploidy for single target conditional on other targets. """

    def __init__(self, mu, covariance, support, data=None, exclude_covar=False):
        self.data = data
        self.support = support
        self.exclude_covar = exclude_covar

        # compute all repeated matrix components for conditional covariance one time upfront
        if not self.exclude_covar:
            self.conditional_mats = [self.get_conditional_mvn(mu, covariance, index, return_matrix_comp=True) for index in xrange(len(mu))]
            logging.info('Finished computing conditional matrix components')

        self.mu = mu
        # keep only diagonal if excluding covariances
        self.covariance = np.diag(np.diagonal(covariance)) if self.exclude_covar else covariance
        self.mu_full = np.concatenate((mu.flatten(), [0]))
        self.inv_covariance_full = np.concatenate((np.concatenate((np.linalg.inv(self.covariance), np.zeros((1,len(self.covariance)))), axis=0),
                                                   np.zeros((len(self.covariance)+1,1))), axis=1)

    def sample(self, copies, intensities, target_index, is_baseline):
        """Given a current set of intensities, and the current ploidy state maintained in this class,
        sample a new ploidy state and intensity for a particular target using the prior distribution as the proposal distribution."""
        intensities_proposed = np.copy(intensities)
        if not is_baseline:
            copies_proposed = np.copy(copies)
            copy_proposed = np.random.choice(self.support, size=1)
            copies_proposed[target_index] = copy_proposed
        else:
            copies_proposed = copies
            copy_proposed = copies[target_index]

        # don't update intensity for last target (for identifiability)
        if target_index != (len(self.mu)):
            if self.exclude_covar:
                mu_bar = self.mu[target_index]
                cov_bar = self.covariance[target_index, target_index]
            else:
                mu_bar, cov_bar = self.get_conditional_mvn(self.mu, self.covariance, target_index,
                                                           intensities[:-1], self.conditional_mats[target_index])

            # sample intensity from conditional normal
            intensity_proposed = np.random.normal(mu_bar, np.sqrt(cov_bar))
            intensities_proposed[target_index] = intensity_proposed

            # calculate log likelihood of proposed jump (J(proposed|current state))
            # Since marginally normal, we just need distance from mean squared
            jump_proposed = -(((intensity_proposed - mu_bar) ** 2) / (2 * cov_bar))
            jump_previous = -(((intensities[target_index] - mu_bar) ** 2) / (2 * cov_bar))
        else:
            intensity_proposed = intensities[target_index]
            jump_proposed = 0
            jump_previous = 0

        # calculate unnormalized log joint likelihood
        joint_proposed = self.log_joint_likelihood(intensities_proposed, copies_proposed)
        joint_previous = self.log_joint_likelihood(intensities, copies)

        log_test_ratio = joint_proposed + jump_previous - joint_previous - jump_proposed

        # sample for new joint state and keep track of acceptance
        # Note python's "or" is equivalent to || not |, saving the RNG and exp
        if log_test_ratio > 0 or np.random.rand() < np.exp(log_test_ratio):
            return copy_proposed, intensity_proposed, 1.0
        else:
            return copies[target_index], intensities[target_index], 0.0

    def log_joint_likelihood(self, intensities, copies, return_neg=False):
        """ Returns unnormalized log likelihood of joint probability given subject data, full set of copy numbers and
            full set of intensities."""
        log_joint = (np.sum(self.data) * -1 * np.log(np.sum(np.multiply(copies, np.exp(intensities)))) +
                     np.sum(np.multiply(self.data, (np.log(copies) + intensities))) +
                     (-0.5 * np.dot(np.dot((intensities - self.mu_full).reshape((1,-1)), self.inv_covariance_full),
                                    (intensities - self.mu_full).reshape((-1,1)))))
        if return_neg:
            log_joint *= -1.
        return log_joint[0][0]

    @staticmethod
    def get_conditional_mvn(mu, cov, index, intensities=None, matrix_comp=None, return_matrix_comp=False):
        """ Returns mu and covariance for conditional normal distribution for single unknown value,
        or covariance matrix components needed to compute conditional mu and covariance.

        mu -- array of len k
        cov -- k x k matrix of covariance values corresponding to mu array
        index -- index of unknown value in intensities
        intensities -- array of len k with known intensity values (including value for unknown index),
                       not required if only returning matrix_comp
        matrix_comp -- dictionary containing pre-computed matrix components,
                       cov_12, cov_22_inv, cov_bar for faster performance
        return_matrix_comp -- will return matrix_comp dictionary, not mu and covariance

        Quadrants are as follows (after translation of desired index to 1,1 position):
        [1,1 [.. ... 1,2 ... ..]]
        [[.. [.. ... ... ... ...]
        [... ... ... ... ... ...]
        [2,1 ... ... 2,2 ... ...]
        [... ... ... ... ... ...]
        [..] ... ... ... ... ..]]

        """

        if matrix_comp is None:
            # split covariance matrix into quadrants and compute each one separately, storing in dict as necessary
            cov_11 = cov[index, index]
            matrix_comp = {'cov_12': np.concatenate((cov[index][:index], cov[index][(index+1):])).reshape((1, -1))}
            cov_21 = np.concatenate((cov[:,index][:index], cov[:, index][(index+1):])).reshape((-1,1))

            cov_22_t = np.concatenate((cov[:index, :index], cov[:index, (index+1):]), axis=1)
            cov_22_b = np.concatenate((cov[(index+1):, :index], cov[(index+1):, (index+1):]), axis=1)
            cov_22 = np.concatenate((cov_22_t, cov_22_b), axis=0)

            matrix_comp['cov_22_inv'] = np.linalg.inv(cov_22)
            matrix_comp['cov_bar'] = cov_11 - np.dot(np.dot(matrix_comp['cov_12'], matrix_comp['cov_22_inv']), cov_21)

        if return_matrix_comp:
            return matrix_comp
        else:
            mu_1 = mu[index]
            mu_2 = np.delete(mu, index)
            a = np.delete(intensities, index)

            mu_bar = mu_1 + np.dot(np.dot(matrix_comp['cov_12'], matrix_comp['cov_22_inv']), (a - mu_2).reshape((-1,1)))
            return mu_bar.flatten()[0], matrix_comp['cov_bar'].flatten()[0]
