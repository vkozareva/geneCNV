import pybedtools



class TargetCollection(object):
    """This class represents an arbitrary list of targets (e.g. DMD exons, TSID bait locations, etc). It
    is useful for organizing/merging targets and can parse a BAM file to produce readcounts across the specified targets"""
    def __init__(self, bedfile_path):
        """Given a list of targets, this should sort them and merge appropriately close intervals """
        self.bed_object = pybedtools.BedTool(bedfile_path)

    def getData(self, bamFileName):
        """Given a BAM file, return a vector with readcounts"""
        pass

    def filter_intervals(self, gene):
        return self.bed_object.filter(lambda x: x[3].split('.')[0] == gene)

    @staticmethod
    def save_intervals(intervals):
        # return [Target(item) for item in intervals]
        return [{'start': intrv.start, 'end': intrv.end, 'chrom': intrv.chrom} for intrv in intervals]
