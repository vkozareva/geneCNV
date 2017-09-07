""" Collection of utilities for generating input to test
that the software works as expected """

import os
import pysam

from numpy.random import normal
from numpy.random import poisson

from cnv.Targets.Target import Target
from cnv.Targets.TargetCollection import TargetCollection


def _validate_file_names(outputPrefix):
    fofn = _make_fofn_name(outputPrefix)
    bams = _make_bam_names(outputPrefix)
    bams.append(fofn)
    for f in bams:
        if os.path.exists(f):
            raise IOError("File named " + f + " was going to be created but already exists.")


def _output_fofn(output_prefix):
    bam_names = _make_bam_names(output_prefix)
    with open(_make_fofn_name(output_prefix), 'w') as f:
        f.write("\n".join(bam_names))


def _make_fofn_name(output_prefix):
    return output_prefix + "files.fofn"


def _make_bam_names(output_prefix, numBams = 10):
    return [output_prefix + "_bam_" + str(x) +'.bam' for x in range(1, (numBams + 1))]

def _write_fake_bam_file(bam_name, header, targets, intensities):
    # Get a dictionary to convert sequence name to index
    name_to_index = {}
    index = 0
    for sq in header["SQ"]:
        name_to_index[sq["SN"]] = index
        index += 1

    with pysam.AlignmentFile(bam_name, "wb", header = header)as bamf:
        rname = 1
        totalCoverage = 30 * len(targets)
        for intensity, t in zip(intensities, targets):
            s = t.start
            e = t.end
            length = max(e - s, 10)
            numreads = poisson(totalCoverage * intensities)[0]
            for replicate in xrange(0, numreads):
                a = pysam.AlignedSegment()
                a.query_name = "Read_" + str(rname)
                a.query_sequence = "A" * length
                a.reference_id = name_to_index[t.chrom]
                bamf.write(a)
                rname += 1




def make_simulated_data(output_prefix, target_file):
    _validate_file_names(output_prefix)
    _output_fofn(output_prefix)

    targets = TargetCollection.load_from_txt_file(target_file)
    header = targets.make_fake_sam_header()

    # Now let's make a fake frequency spectrum
    intensities = abs(normal(size=len(targets)))
    intensities = intensities / sum(intensities)
    # Write a bunch of bam files
    for bam in _make_bam_names(output_prefix):
        _write_fake_bam_file(bam, header, targets, intensities)



