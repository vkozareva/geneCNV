import os
from Target import Target


class TargetCollection(list):
    """This class represents an arbitrary list of targets (e.g. DMD exons, TSID bait locations, etc). It
    is useful for organizing/merging targets and can parse a BAM file to produce readcounts across the specified targets"""
    def __init__(self, targets = None):
        """
        Initialize a Target Collection with a list of targets, each of which should be
        :param targets: a Collection
        """
        super(TargetCollection, self).__init__()
        if targets:
            for t in targets:
                if not isinstance(t, Target):
                    raise TypeError("Tried to create target collection with non target types.")
                self.append(t)
        self.sort()

    def make_fake_sam_header(self):
        """
        Makes a header to create a SAM/BAM file compatible with the regions defined
        in this program.  Used for simulating BAM files.

        :return: The dictionary with header information
        """
        header = {"HD": {"VN": "1.0"}}
        # Now we'll make fake sequences by iterating through each
        # Target and finding the maximum of the range
        seqs = []
        cur_seq = self[0].chrom
        cur_end = 0
        for t in self:
            if t.chrom == cur_seq:
                cur_end = max(cur_end, t.end)
            else:
                seqs.append({"LN": cur_end+1,
                             "SN": cur_seq})
        seqs.append({"LN": cur_end + 1,
                     "SN": cur_seq})
        header["SQ"] = seqs
        return header

    def write_to_txt_file(self, fname):
        if os.path.exists(fname):
            raise IOError("Output filename already exists")
        with open(fname, 'w') as f:
            for t in self:
                f.write(str(t) + "\n")


    @classmethod
    def load_from_txt_file(cls, fname, query_interval=None):
        """
        Loads a plain text bed file into a list of Interval objects,
        only including those in the query interval if specified.

        :param fname: text file with bed format of chrom, start, end, id
        :param query_interval: Query Interval, will return all if none
        :return: A TargetCollection from the bed file that overlap the query interval
        """
        if not os.path.exists(fname):
            raise IOError("Text file " + fname + " does not exist")
        result = cls()
        with open(fname) as fh:
            for line in fh:
                if line[0] != "#":
                    interval = Target.create_from_BED_string(line)
                    if query_interval and query_interval.overlaps(interval):
                        result.append(interval)
                    elif not query_interval:
                        result.append(interval)
        result.sort()
        return result
