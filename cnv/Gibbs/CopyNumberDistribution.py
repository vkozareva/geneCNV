from cnv.Targets.TargetCollection import TargetCollection

class CopyNumberDistribution(object):
    """A class which describes the distribution of ploidy across the different exons"""
    def __init__(self, targets, data):
        """Initializes the data to hold past and current samples of ploidy levels for all values in targets"""
        isinstance(targets, TargetCollection)

    def sample(self, intensities):
        """Given a current set of intensities, and the current ploidy state maintained in this class,
        sample a new ploidy state """
        pass

    def output(self, out_file_name):
        """Report on the posterior distribution obtained by the sampling procedure"""
        pass
