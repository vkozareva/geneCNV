import os
import shutil
import tempfile
import unittest
import pysam

from cnv.utilities import SimulateData
from cnv.Targets.Target import Target
from cnv.Targets.TargetCollection import TargetCollection
from cnv.coverage_matrix import *

def _makeFakeTargets():
    t1 = Target("1", 10, 200, "exon1")
    t2 = Target("1", 1000, 1100, "exon2")
    t3 = Target("2", 10, 200, "exon3")
    return TargetCollection([t1, t2, t3])

class SimulatorTest(unittest.TestCase):
    def test_simulator(self):
        targs = _makeFakeTargets()
        outdir = tempfile.mkdtemp() + "/"
        bedFile = tempfile.mktemp()
        targs.write_to_txt_file(bedFile)
        SimulateData.make_simulated_data(outdir, bedFile)

        expected_names = SimulateData._make_bam_names(outdir)
        for n in expected_names:
            self.assertTrue(os.path.exists(n))
            os.remove(n)

        expected_fofn = SimulateData._make_fofn_name(outdir)
        self.assertTrue(os.path.exists(expected_fofn))


        # Verify we can read those files
        cm = CoverageMatrix()
        cm.create_coverage_matrix(expected_fofn, targs)
        os.remove(expected_fofn)
