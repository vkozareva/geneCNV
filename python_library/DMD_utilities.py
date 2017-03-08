from genepeeks.common import utilities as util


def get_DMD_exons_merged(min_dist=200, only_primary_transcript=True, only_coding_regions=False, print_merged_exons=True):
    """ Get DMD exons from primary DMD transcript in mongo, and merge all exons within min_dist of each other """
    DMD_ensembl = util.Mongo.get_collection_data('gene', wanted_db='prod', query={'_id': 'DMD'}, find_one=True,
                                                 single_field='ensembl', required=True, print_stats=False)
    if only_primary_transcript:
        DMD_primary_transcript_d = util.get_nested_value(DMD_ensembl, ('is_primary', 'transcripts', 'is_primary'), required=True)
        DMD_exons = util.get_transcript_coding_regions(DMD_primary_transcript_d) if only_coding_regions else DMD_primary_transcript_d['exons']
    else:
        # If not only_primary_transcript, combine exon info from all transcripts. Optionally restrict to just coding_regions
        DMD_exons = []
        for transcript_d in DMD_ensembl[0]['transcripts']:
            DMD_exons += util.get_transcript_coding_regions(transcript_d) if only_coding_regions else transcript_d['exons']
        DMD_exons = util.merge_intervals(DMD_exons)

    # Merge all exons within min_dist of each other and use exon_labels to keep track of which exons were merged
    DMD_exons_merged = util.merge_intervals(DMD_exons, min_dist=min_dist, include_index=True)
    exon_labels = ['Ex' + exon['index'] for exon in DMD_exons_merged]
    if print_merged_exons:
        print 'There are {} merged exons, with the following exons being merged: {}'.format(
            len(DMD_exons_merged), [exon['index'] for exon in DMD_exons_merged if not exon['index'].isdigit()])
    return DMD_exons_merged, exon_labels
