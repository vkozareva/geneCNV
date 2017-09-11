import unittest
from cnv.Targets.Target import Target
from cnv.Targets.TargetCollection import TargetCollection


class TargetTests(unittest.TestCase):
    def test_targetSort(self):

        rs = [Target("C", 50, 100),
              Target("A", 75, 100),
              Target("A", 75, 150),
              Target("A", 80, 90)]
        old = list(rs)
        rs.sort()
        order = [2, 1, 3, 0]
        for idx, o in enumerate(order):
            self.assertEquals(rs[idx].start, old[o].start, "Start value was off")
            self.assertEquals(rs[idx].end, old[o].end, "End value was off")
            self.assertEquals(rs[idx].chrom, old[o].chrom, "Start value was off")

    def test_can_enumerate(self):
        rs = [Target("C", 50, 100),
              Target("A", 750, 1500),
              Target("A", 80, 90)]
        nt = TargetCollection(rs)
        q = nt[1:3]
        v= list(enumerate(q))
        self.assertTrue(len(v) == 2)
    def test_targetcollection_merge(self):
        rs = [Target("C", 50, 100),
              Target("A", 75, 100, "first"),
              Target("A", 75, 150),
              Target("A", 80, 90)]
        # overlapping intervals not allowed
        self.assertRaises(lambda : TargetCollection(rs))

        # Now verify non-overlapping intervals on different chromosomes are merged.
        nl = rs[:2]
        nl.append(Target("A", 300, 400, "second"))
        nt = TargetCollection(nl)
        merged = nt[0]
        self.assertTrue(len(nt) == 2)
        self.assertTrue(merged.start == 75)
        self.assertTrue(merged.end == 400)
        self.assertTrue(merged.chrom == "A")
        self.assertTrue(merged.label.startswith("Merged"))


    def test_targetOverlaps(self):
        a = Target("A", 150, 200)
        b = Target("B", 150, 200)
        self.assertFalse(a.overlaps(b))
        c = Target("A", 200, 210)
        self.assertFalse(a.overlaps(c))
        c.start = 199
        self.assertTrue(a.overlaps(c))
        i1 = Target("1", 10, 20, "")
        i2 = Target("1", 19, 30, "")
        i3 = Target("2", 200, 5000, "")
        i4 = Target("1", 5000, 10000)
        i5 = Target("1", 20, 25)
        i6 = Target("1", 12, 17)
        self.assertTrue(i1.overlaps(i2))
        self.assertTrue(i2.overlaps(i1))
        self.assertFalse(i3.overlaps(i1))
        self.assertFalse(i4.overlaps(i1))
        self.assertFalse(i5.overlaps(i1))
        self.assertTrue(i6.overlaps(i1))

if __name__ == '__main__':
    unittest.main()

