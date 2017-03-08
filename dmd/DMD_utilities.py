from genepeeks.common import utilities as util


def get_merged_exons(min_dist=200, only_primary_transcript=True, only_coding_regions=False, print_merged_exons=True, wanted_gene='DMD'):
    """ Get DMD exons from primary DMD transcript in mongo, and merge all exons within min_dist of each other """
    ensembl_data = util.Mongo.get_collection_data('gene', wanted_db='prod', query={'_id': wanted_gene}, find_one=True,
                                                  single_field='ensembl', required=True, print_stats=False)
    if only_primary_transcript:
        primary_transcript_d = util.get_nested_value(ensembl_data, ('is_primary', 'transcripts', 'is_primary'), required=True)
        exons = util.get_transcript_coding_regions(primary_transcript_d) if only_coding_regions else primary_transcript_d['exons']
    else:
        # If not only_primary_transcript, combine exon info from all transcripts. Optionally restrict to just coding_regions
        exons = []
        for transcript_d in ensembl_data[0]['transcripts']:
            exons += util.get_transcript_coding_regions(transcript_d) if only_coding_regions else transcript_d['exons']
        exons = util.merge_intervals(exons)

    # Merge all exons within min_dist of each other and use exon_labels to keep track of which exons were merged
    merged_exons = util.merge_intervals(exons, min_dist=min_dist, include_index=True)
    exon_labels = ['Ex' + exon['index'] for exon in merged_exons]
    if print_merged_exons:
        print 'There are {} merged exons, with the following exons being merged: {}'.format(
            len(merged_exons), [exon['index'] for exon in merged_exons if not exon['index'].isdigit()])
    return merged_exons, exon_labels
