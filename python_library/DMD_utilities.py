from genepeeks.common import utilities as util


def get_DMD_exons_merged(min_dist=200, print_merged_exons=True):
    """ Get DMD exons from primary DMD transcript in mongo, and merge all exons within min_dist of each other """
    DMD_ensembl = util.Mongo.get_collection_data('gene', wanted_db='prod', query={'_id': 'DMD'}, find_one=True,
                                                 single_field='ensembl', required=True, print_stats=False)
    DMD_primary_transcript_d = util.get_nested_value(DMD_ensembl, ('is_primary', 'transcripts', 'is_primary'), required=True)
    DMD_exons = DMD_primary_transcript_d['exons']

    # Add coding start and end to first and last exons
    DMD_coding_region = DMD_primary_transcript_d['translation']
    DMD_exons[0]['coding_start'] = DMD_coding_region['start']
    DMD_exons[-1]['coding_end'] = DMD_coding_region['end']

    # Merge all exons within min_dist of each other and use exon_labels to keep track of which exons were merged
    DMD_exons_merged = util.merge_intervals(DMD_exons, min_dist=min_dist, include_index=True)
    exon_labels = ['Ex' + exon['index'] for exon in DMD_exons_merged]
    if print_merged_exons:
        print 'The following exons were merged: {}'.format([exon['index'] for exon in DMD_exons_merged if not exon['index'].isdigit()])
    return DMD_exons_merged, exon_labels


def get_exon_num(read_start, read_end, exon_list, skipped_counts):
    """ Determine which exon a read overlaps with, if any """
    exon_num = None
    for i, exon in enumerate(exon_list):
        if exon['end'] < read_start:
            continue
        else:
            if exon['start'] <= read_end:
                # Check if the read pair falls in two exons, and skip if it does.
                if i != (len(exon_list) - 1) and exon_list[i + 1]['start'] <= read_end:
                    util.add_to_dict(skipped_counts, 'in_two_exons')
                else:
                    exon_num = i
            else:
                util.add_to_dict(skipped_counts, 'outside_of_exon')
            break
    return exon_num
