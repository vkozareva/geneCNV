

class Target(object):
    """A simple class which holds a chromosome region."""
    __slots__ = ('chromosome', 'start', 'end')

    def __init__(self, chromosome, start, end):
        """Create a new target region with start (inclusive) and end (exclusive)"""
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def overlaps(self, other):
        assert isinstance(other, Target)
        return self.chromosome == other.chromosome and self.start < other.end and self.end > other.start

    def __cmp__(self, other):
        """Orders by name, start, then largest length first"""
        assert isinstance(other, Target)
        chrom_comp = cmp(self.chromosome, other.chromosome)
        if chrom_comp == 0:
            start_cmp = cmp(self.start, other.start)
            if start_cmp == 0:
                return cmp(other.end, self.end) # Note, longest interval goes first
            else:
                return start_cmp
        else:
            return chrom_comp

    def __repr__(self):
        return '-'.join([self.chromosome, str(self.start), str(self.end)])

    # Python 3 boilerplate
    def __lt__(self, other):
        return self.__cmp__(other) < 0

    def __gt__(self, other):
        return self.__cmp__(other) > 0

    def __eq__(self, other):
        return self.__cmp__(other.obj) == 0

    def __le__(self, other):
        return self.__cmp__(other.obj) <= 0

    def __ge__(self, other):
        return self.__cmp__(other.obj) >= 0

    def __ne__(self, other):
        return self.__cmp__(other.obj) != 0