import numpy as np

class IntensitiesDistribution(object):
    """A class which stores our current parameters and contains methods to update them during Gibbs Steps"""
    def __init__(self, parametersFile, scale):
        self.scale = scale

    def sample(self, X_priors, data=None):
        """Given a current ploidy state, and an internally maintained vector of current intensities, sample a
        new intensity vector.  This will be a no-op until we can better define the model."""
        return (np.random.dirichlet(X_priors * self.scale + data) if data is not None
                else np.random.dirichlet(X_priors * self.scale))
