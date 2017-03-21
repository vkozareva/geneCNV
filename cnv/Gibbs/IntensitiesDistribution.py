class IntensitiesDistribution(object):
    """A class which stores our current parameters and contains methods to update them during Gibbs Steps"""
    def __init__(self, parameterFileName=None, intensities=None):
        self.parameterFileName = parameterFileName
        self.intensities = intensities

    def sample(self, copyNumberModel, data=None):
        """Given a current ploidy state, and an internally maintained vector of current intensities, sample a
        new intensity vector.  This will be a no-op until we can better define the model."""
        pass
