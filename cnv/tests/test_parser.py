import unittest
import pysam

from test_resources import *


class ParserTest(unittest.TestCase):
    def test_dataFileExists(self):
        bamFile = EXAMPLE_BAM_PATH
        self.assertTrue(os.path.exists(bamFile))

    def test_pySamWorks(self):
        bamfile = pysam.AlignmentFile(EXAMPLE_BAM_PATH, "rb")
        count = 0
        # Get coverage data for each sample within each exon
        for read in bamfile.fetch('X', start=31137344, end=33357726):
            count += 1
        # Confirm it matches result from:
        # samtools view example.bam X:31137344-33357726 | wc -l
        self.assertEquals(count, 14518)

if __name__ == '__main__':
    unittest.main()
