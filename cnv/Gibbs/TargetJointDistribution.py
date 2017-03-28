import numpy as np
from IntensitiesDistribution import IntensitiesDistribution
from CopyNumberDistribution import CopyNumberDistribution

class TargetJointDistribution(object):
    """Describe the distribution of ploidy across a set of targets.
    Random data is substituted if optional subject data is omitted.
    Includes methods for Gibbs Sampling of the posterior likelihood."""
    def __init__(self, targets, mu, covariance, data=None, support=None, include_covar=False):
        """Initialize the data to hold past and current samples of ploidy levels for all values in targets."""
        # isinstance(targets, TargetCollection)
        self.data = data
        if support is not None:
            self.support = support
        else:
            self.support = [1, 2, 3]
        # self.copies = copies
        # self.intensities = intensities
        self.include_covar = include_covar
        if self.include_covar:
            self.cond_inverses = np.zeros((len(mu), len(mu)-1, len(mu)-1))
            for index in xrange(len(mu)):
                cov_22_t = np.concatenate((covariance[:index, :index], covariance[:index, (index+1):]), axis=1)
                cov_22_b = np.concatenate((covariance[(index+1):, :index], covariance[(index+1):, (index+1):]), axis=1)
                cov_22 = np.concatenate((cov_22_t, cov_22_b), axis=0)
                cov_22_inv = np.linalg.inv(cov_22)
                self.cond_inverses[index] = cov_22_inv

        self.mu = mu
        self.covariance = covariance
        self.mu_full = np.concatenate((mu.flatten(), [0]))
        self.inv_covariance_full = np.concatenate((np.concatenate((np.linalg.inv(covariance), np.zeros((1,len(covariance)))), axis=0),
                                                np.zeros((len(covariance)+1,1))), axis=1)

    def sample(self, copies, intensities, target_index):
        """Given a current set of intensities, and the current ploidy state maintained in this class,
        sample a new ploidy state and intensity."""
        intensities_proposed = np.copy(intensities)
        copies_proposed = np.copy(copies)

        copy_proposed = CopyNumberDistribution(1, self.support).sample()
        copies_proposed[target_index] = copy_proposed

        # don't update intensity for last target (for identifiability)
        if target_index != (len(self.mu)):
            if self.include_covar:
                mu_bar = self.mu_full[target_index]
                cov_bar = self.covariance[target_index, target_index]
            else:
                mu_bar, cov_bar = self.get_conditional_mvn(self.mu, self.covariance, target_index,
                                                           intensities[:-1], self.cond_inverses[target_index])

            intensity_proposed = np.random.normal(mu_bar, np.sqrt(cov_bar))
            intensities_proposed[target_index] = intensity_proposed

            # calculate
            jump_proposed = -(intensity_proposed - mu_bar) ** 2 / (2 * cov_bar)
            jump_previous = -(intensities[target_index] - mu_bar) ** 2 / (2 * cov_bar)
        else:
            jump_proposed = 0
            jump_previous = 0

        joint_proposed = self.log_joint_likelihood(copies_proposed, intensities_proposed)
        joint_previous = self.log_joint_likelihood(copies, intensities)

        log_test_ratio = joint_proposed + jump_previous - joint_previous - jump_proposed

        # sample for new joint state and keep track of acceptance
        jump_chance = min(1.0, np.exp(log_test_ratio))
        if np.random.rand() < jump_chance:
            return copy_proposed, intensity_proposed, 1
        else:
            return copies[target_index], intensities[target_index], 0

        # cnv_probs = np.zeros((len(intensities.intensities), len(self.support)))
        # for exon in range(len(intensities.intensities)):
        #     cnv_log_probs = np.zeros(len(self.support))
        #     for value in self.support:
        #         self.copies[exon] = value
        #         # get new normed probabilities given test value and priors for exon intensities
        #         cnv_log_probs[value - 1] = self.log_likelihood(intensities.intensities)
        #     cnv_log_probs = cnv_log_probs - np.max(cnv_log_probs)
        #     cnv_probs[exon] = np.exp(cnv_log_probs)
        #     cnv_probs[exon] = cnv_probs[exon] / np.sum(cnv_probs[exon])
        #     new_copies = np.random.choice(self.support, p = cnv_probs[exon])
        #     self.copies[exon] = new_copies

        # return self.log_likelihood(intensities.intensities), cnv_probs

    def log_joint_likelihood(self, copies, intensities):
        log_joint = (np.sum(self.data) * -1 * np.log(np.sum(np.multiply(copies, np.exp(intensities)))) +
                     np.sum(np.multiply(self.data, (np.log(copies) + intensities))) +
                     (-0.5 * np.dot(np.dot((intensities - self.mu_full).reshape((1,-1)), self.inv_covariance_full),
                                    (intensities - self.mu_full).reshape((-1,1)))))
        return log_joint

    @staticmethod
    def get_conditional_mvn(mu, cov, index, intensities, cov_22_inv=None):
        """Returns mu and cov for conditional normal distribution for single unknown value
           Pass in previously inverted 2,2 component for much faster performance"""
        mu_1 = mu[index]
        mu_2 = np.delete(mu, index)
        a = np.delete(full_x, index)
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


    # def log_likelihood(self, X_probs):
    #     normed_probs = self.copies * X_probs
    #     normed_probs /= np.sum(normed_probs)
    #     return np.sum(np.log(normed_probs) * self.data)