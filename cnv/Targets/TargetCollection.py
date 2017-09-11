import os
from Target import Target

DEFAULT_MERGE_DISTANCE = 629


class TargetCollection(object):
    """This class represents an arbitrary list of targets (e.g. DMD exons, TSID bait locations, etc). It
    is useful for organizing/merging targets. """
    def __init__(self, targets = None, min_merge_dist = DEFAULT_MERGE_DISTANCE):
        """
        Initialize a TargetCollection with a list of targets, each of which should be a target
        :param targets: a Collection
        :param min_merge_dist: The minimum distance that nearby intervals will be merged with
        """
        self._collection = list()
        self._sorted = False
        self._merged = False
        self.min_dist = min_merge_dist

        if targets:
            for t in targets:
                self.append(t)
        if self._collection: # Empty sequences are false
            self.merge_nearby_intervals(min_merge_dist)


    def merge_nearby_intervals(self, min_dist = DEFAULT_MERGE_DISTANCE):
        """
        For the intervals in this list, confirm that non currently overlap, and merge any
        that are within the min_dist specified on the same chromosome.

        :param min_dist: Intervals <= this amount will be merged.
        :return:
        """
        if self._merged:
            return
        self.sort()
        self.min_dist = min_dist

        # Confirm no overlaps in the intervals and get distances between them
        cur_chrom = self._collection[0].chrom
        cur_end = self._collection[0].end
        difs = []
        DIF_CHROM_FLAG = -1
        for t in self._collection[1:]:
            if t.chrom == cur_chrom:
                if t.start <= cur_end:
                    raise IOError("TargetCollection contains overlapping intervals")
                difs.append(t.start - cur_end)
            else:
                difs.append(DIF_CHROM_FLAG)
            cur_chrom = t.chrom
            cur_end = t.end

        # Now merge any intervals that were smaller than our min interval
        # assume this is rare so modifying list in place instead of copying.
        difs.reverse()
        cur_top = len(self._collection) - 1
        for space in difs:
            if space != DIF_CHROM_FLAG and space <= min_dist:
                top = self._collection[cur_top]
                bottom = self._collection[cur_top - 1]
                new_target = bottom.Target(top)
                self._collection.pop(cur_top)
                self._collection[cur_top - 1] = new_target
            cur_top -= 1
        self._merged = True


    def append(self, item):
        if not isinstance(item, Target):
            raise TypeError("TargetCollection requires Target objects")
        self._sorted = False
        self._merged = False
        self._collection.append(item)

    def sort(self):
        if not self._sorted:
            self._collection.sort()
            self._sorted = True

    def t_index(self, item):
        if not isinstance(item, Target):
            raise TypeError('Object to find must be Target')
        return self._collection.index(item)

    def make_fake_sam_header(self):
        """
        Makes a header to create a SAM/BAM file compatible with the regions defined
        in this program.  Used for simulating BAM files.

        :return: The dictionary with header information
        """
        header = {"HD": {"VN": "1.0"}}
        self.sort()
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
                cur_end = t.end
                cur_seq = t.chrom
        seqs.append({"LN": cur_end + 1,
                     "SN": cur_seq})
        header["SQ"] = seqs
        return header

    def write_to_txt_file(self, fname):
        if os.path.exists(fname):
            raise IOError("Output filename already exists")
        self.sort()
        with open(fname, 'w') as f:
            for t in self:
                f.write(str(t) + "\n")


    @classmethod
    def load_from_txt_file(cls, fname, query_interval=None, min_merge_dist=DEFAULT_MERGE_DISTANCE):
        """
        Loads a plain text bed file into a list of Interval objects,
        only including those in the query interval if specified.

        :param fname: text file with bed format of chrom, start, end, id
        :param query_interval: Query Interval, will return all if none
        :param min_merge_dist: Minimum distance to merge intervals by
        :return: A TargetCollection from the bed file that overlap the query interval
        """
        if not os.path.exists(fname):
            raise IOError("Text file " + fname + " does not exist")
        result = cls(min_merge_dist = min_merge_dist)
        with open(fname) as fh:
            for line in fh:
                if line[0] != "#":
                    interval = Target.create_from_BED_string(line)
                    if query_interval and query_interval.overlaps(interval):
                        result.append(interval)
                    elif not query_interval:
                        result.append(interval)
        result.merge_nearby_intervals()
        return result

    ## List method implementations, we override this to guarantee we maintain
    # the sorted and merged state
    def __len__(self):
        return len(self._collection)

    def __getitem__(self, i):
        if isinstance(i, slice):
            newtargets = self._collection[i]
            # Get the start, stop, and step from the slice
            return TargetCollection(newtargets, min_merge_dist=self.min_dist)
        return self._collection[i]

    def __delitem__(self, i):
        del self._collection[i]

    def __setitem__(self, key, value):
        self._collection[key] = value
        self._sorted = False
        self._merged = False

    def __iter__(self):
        return self._collection.__iter__()