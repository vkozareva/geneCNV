import numpy as np

class IntensitiesDistribution(object):
    """A class which stores our current parameters and contains methods to update them during Gibbs Steps"""
    def __init__(self, mu, covariance):
        self.mu = mu
        self.covariance = covariance

    def sample(self):
        """Given a current ploidy state, and an internally maintained vector of current intensities, sample a
        new intensity vector.  This will be a no-op until we can better define the model."""
        intensities = np.concatenate((np.random.multivariate_normal(self.mu.flatten(), self.covariance), [0]))
        return intensities
