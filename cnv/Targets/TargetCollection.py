from Target import Target

class TargetCollection(object):
    """This class represents an arbitrary list of targets (e.g. DMD exons, TSID bait locations, etc). It
    is useful for organizing/merging targets and can parse a BAM file to produce readcounts across the specified targets"""
    def __init__(self, targets):
        """Given a list of targets, this should sort them and merge appropriately close intervals """
        pass
    def getData(self, bamFileName):
        """Given a BAM file, return a vector with readcounts"""
        pass