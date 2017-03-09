import unittest
from cnv.Targets.Target import Target

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
            self.assertEquals(rs[idx].chromosome, old[o].chromosome, "Start value was off")

    def test_targetOverlaps(self):
        a = Target("A", 150, 200)
        b = Target("B", 150, 200)
        self.assertFalse(a.overlaps(b))
        c = Target("A", 200, 210)
        self.assertFalse(a.overlaps(c))
        c.start = 199
        self.assertTrue(a.overlaps(c))
