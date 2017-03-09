"""The main Gibbs Sampling Model """
from IntensitiesDistribution import IntensitiesDistribution
from CopyNumberDistribution import CopyNumberDistribution
from cnv.Targets.TargetCollection import TargetCollection

class PloidyModel(object):
    """This is the full statistical model and class that runs the Gibbs Sampling. It is responsible for taking a
    parameter set, a BAM file, and a TargetCollection and running a Gibbs sampler to determine the
     posterior probability of different ploidy states."""

    def __init__(self, bamFileName, targets, parameterFileName):
        """Initialize the data model with it's input arguments, we should load the parameters and
        calculate the coverage at each interval in the BAM file"""
        isinstance(targets, TargetCollection)
        # Fake code below
        data = TargetCollection.getData(bamFileName)
        self.intensities = IntensitiesDistribution.IntensitiesDistribution(parameterFileName)
        self.ploidy = CopyNumberDistribution(targets, data)

    def RunGibbsSampler(self, numIterations):
        """ Run a Gibbs sampler for several iterations """
        for i in xrange(0, numIterations):
            self.intensities.sample(self.ploidy)
            self.ploidy.sample(self.intensities)

    def Output(self, fileName):
        """Output a results file and a PDF of the posterior probabilities plot."""
        self.ploidy.output(fileName)
