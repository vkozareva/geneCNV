import pickle

class HLN_Parameters(object):
    """A simple class for serializing/deserializing the parameters needed for subject testing.
    This currently includes the target intervals, and the model hyperparameters."""
    def __init__(self, targets, mu, covariance):
        self.targets = targets
        self.mu = mu
        self.covariance = covariance

    def dump(self, pickleFile):
        pickle.dump(self.targets, pickleFile)
        pickle.dump(self.mu, pickleFile)
        pickle.dump(self.covariance, pickleFile)

    def load(self, pickleFile):
        self.targets = pickle.load(pickleFile)
        self.mu = pickle.load(pickleFile)
        self.covariance = pickle.load(pickleFile)
