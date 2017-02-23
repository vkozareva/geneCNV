from genepeeks.common import utilities as util


def get_DMD_exons_merged(min_dist=200):
    DMD_ensembl = util.Mongo.get_collection_data('gene', wanted_db='prod', query={'_id': 'DMD'}, find_one=True, single_field='ensembl')
    DMD_exons = util.get_nested_value(DMD_ensembl, ('is_primary', 'transcripts', 'is_primary', 'exons'))
    DMD_exons_merged = util.merge_intervals(DMD_exons, min_dist=min_dist, include_index=True)
    exon_labels = ['Ex' + exon['index'] for exon in DMD_exons_merged]
    # print len(DMD_exons_merged), [exon['index'] for exon in DMD_exons_merged if not exon['index'].isdigit()]
    return DMD_exons_merged, exon_labels


def get_exon_num(read_start, read_end, DMD_exons_merged, skipped_counts):
    """ Determine which exon a read overlaps with, if any """
    exon_num = None
    for i, exon in enumerate(DMD_exons_merged):
        if exon['end'] <= read_start:
            continue
        else:
            if exon['start'] < read_end:
                exon_num = i
                # If not already in the last exon, check if the read pair falls in the next exon as well.
                if exon_num != (len(DMD_exons_merged) - 1):
                    next_exon = DMD_exons_merged[i + 1]
                    if next_exon['start'] < read_end:
                        util.add_to_dict(skipped_counts, 'in_two_exons')
                        # Skip the read pair if it falls in both exons
                        exon_num = None
            else:
                util.add_to_dict(skipped_counts, 'outside_of_exon')
            break
    return exon_num
