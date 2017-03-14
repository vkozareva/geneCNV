from genepeeks.common import utilities as util
from Targets.TargetCollection import TargetCollection
import os


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


def get_interval_labels(combined_intervals_merged, wanted_gene):
    """ Get labels for each interval. If interval overlaps with an exon in the primary transcript,
    use primary transcript exon number. Otherwise use target number """

    ensembl_data = util.Mongo.get_collection_data('gene', wanted_db='prod', query={'_id': wanted_gene}, find_one=True,
                                                  single_field='ensembl', required=True, print_stats=False)

    primary_exons = util.get_nested_value(ensembl_data, ('is_primary', 'transcripts', 'is_primary', 'exons'), required=True)
    interval_labels = []
    for i, intrv in enumerate(combined_intervals_merged):
        # Check if interval overlaps with the primary transcript, and use that exon number if it does
        primary_exon_indexes = util.in_interval((intrv['start'], intrv['end']), primary_exons)
        if primary_exon_indexes:
            label = 'Ex{}'.format(','.join([str(index + 1) for index in primary_exon_indexes]))
        else:
            label = 'Target{}'.format(i + 1)
        interval_labels.append(label)
    return interval_labels


def combine_panel_intervals(wanted_gene='DMD', min_dist=629):
    """ Get interval union of intervals from TruSight One and TruSight Inherited Disease """
    bedfile_paths = {
        'TSID': {'file': os.path.join('..', 'inputs', 'TruSight_Inherited_Disease_Manifest_A.bed')},
        'TSO': {'file': os.path.join('..', 'inputs', 'TruSight-One-BED-May-2014.txt')}
    }

    intervals_list = []
    # Convert each bed file to a list of intervals for the desired gene
    for name, intrv_info in bedfile_paths.items():
        instance = TargetCollection(intrv_info['file'])
        filtered_intervals = instance.filter_intervals(wanted_gene)
        intervals_list.append(instance.save_intervals(filtered_intervals))

    # Get union of the intervals and then merge
    combined_intervals = util.interval_union(*intervals_list)
    combined_intervals_merged = util.merge_intervals(combined_intervals, min_dist=min_dist)

    interval_labels = get_interval_labels(combined_intervals_merged, wanted_gene)

    return combined_intervals_merged, interval_labels
