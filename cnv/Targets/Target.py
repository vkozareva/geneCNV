import sys

class Target(object):
    """A simple class which holds a genomic region."""
    __slots__ = ('chrom', 'start', 'end', 'label', 'name')

    def __init__(self, chrom, s=0, e=sys.maxint, label="", name=None):
        """
        Create a new target region with start (inclusive) and end (exclusive)

        :param chrom: string The chromosome, e.g. 'chrX'
        :param s: int start
        :param e: int end
        :param label: The id of the target
        :param name: string An additional label, usually the gene name
        """

        self.chrom = chrom
        self.start = s
        self.end = e
        self.label = label
        self.name = name

    def overlaps(self, other):
        assert isinstance(other, Target)
        # assert self.chrom == other.chrom
        # return min(self.end, other.end) - max(self.start, other.start)
        return self.chrom == other.chrom and self.start < other.end and self.end > other.start

    def merge(self, other):
        """ Merges in another Target object, concatenating the names at least once,
        and if more than one merge simply renaming the new interval as Merged-START:END
        :param other the Target to merge in.
        :returns A new target with the full extent of both intervals
        """
        assert isinstance(other, Target)
        assert other.chrom == self.chrom
        s = min(self.start, other.start)
        e = max(self.end, other.end)
        mergedString = "MergedBaseline-" if "Baseline" in self.label or "Baseline" in other.label else "Merged-"
        if self.label.count(mergedString) > 0 or other.label.count(mergedString) > 0:
            new_name = mergedString + str(s) + ":" + str(e)
        else:
            new_name = mergedString + self.label + ":" + other.label
        return Target(self.chrom, s, e, new_name)

    def __cmp__(self, other):
        """Orders by chrom, start, then largest length first
        Baseline targets will always be second"""
        assert isinstance(other, Target)
        # for easier chrom comparison
        convert = lambda chrom: int(chrom) if chrom.isdigit() else chrom

        base_comp = cmp('Baseline' in self.label, 'Baseline' in other.label)
        if base_comp == 0:
            chrom_comp = cmp(convert(self.chrom), convert(other.chrom))
            if chrom_comp == 0:
                start_cmp = cmp(self.start, other.start)
                if start_cmp == 0:
                    return cmp(other.end, self.end)     # Note, longest interval goes first

                return start_cmp

            return chrom_comp

        return base_comp

    def __str__(self):
        return "\t".join([self.chrom, str(self.start),
                          str(self.end), self.label])
    def __repr__(self):
        return '{}:{}-{}'.format(self.chrom, self.start, self.end)

    # Python 3 boilerplate
    # In python 3 we can't simply implement the __cmp__ method but must implement
    # all of these rich comparisons, added for forward compatibility
    def __lt__(self, other):
        return self.__cmp__(other) < 0

    def __gt__(self, other):
        return self.__cmp__(other) > 0

    def __eq__(self, other):
        return self.__cmp__(other) == 0

    def __le__(self, other):
        return self.__cmp__(other) <= 0

    def __ge__(self, other):
        return self.__cmp__(other) >= 0

    def __ne__(self, other):
        return self.__cmp__(other) != 0

    @classmethod
    def create_from_BED_string(cls, line):
        sp = line.strip().split("\t")
        if len(sp) == 5:
            return cls(sp[0], int(sp[1]), int(sp[2]), sp[3], sp[4])
        elif len(sp) == 4:
            return cls(sp[0], int(sp[1]), int(sp[2]), sp[3])
        elif len(sp) == 3:
            return cls(sp[0], int(sp[1]), int(sp[2]))
        else:
            raise RuntimeError("Text file line did not contain three fields of CHROM, START,END:\n " +
                               "Error on line: " + line)
