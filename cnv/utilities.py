import os
from Targets.TargetCollection import TargetCollection
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


def add_interval_labels(combined_intervals_merged, primary_exons):
    """ Get labels for each interval. If interval overlaps with an exon in the primary transcript,
    use primary transcript exon number. Otherwise use target number """

    for i, intrv in enumerate(combined_intervals_merged):
        # Check if interval overlaps with the primary transcript, and use that exon number if it does
        primary_exon_indexes = util.in_interval((intrv['start'], intrv['end']), primary_exons)
        if primary_exon_indexes:
            label = 'Ex{}'.format('-'.join([str(index + 1) for index in primary_exon_indexes]))
        else:
            label = 'Target{}'.format(i + 1)
        intrv['label'] = label


def combine_panel_intervals(wanted_gene='DMD', min_dist=629):
    """ Get interval union of intervals from TruSight One and TruSight Inherited Disease """
    input_dir = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'inputs')
    bedfile_paths = {
        'TSID': {'file': os.path.join(input_dir, 'TruSight_Inherited_Disease_Manifest_A.bed')},
        'TSO': {'file': os.path.join(input_dir, 'TruSight-One-BED-May-2014.txt')},
        'primary': {'file': os.path.join(input_dir, 'primary_{}_exons.bed'.format(wanted_gene))},
    }

    # For now, commenting out the code that takes the primary exons from Mongo and instead taking them from a saved bed file
    # ensembl_data = util.Mongo.get_collection_data('gene', wanted_db='prod', query={'_id': wanted_gene}, find_one=True,
    #                                               single_field='ensembl', required=True, print_stats=False)
    # primary_exons = util.get_nested_value(ensembl_data, ('is_primary', 'transcripts', 'is_primary', 'exons'), required=True)

    intervals_list = []
    # Convert each bed file to a list of intervals for the desired gene
    for name, intrv_info in bedfile_paths.items():
        if not os.path.exists(intrv_info['file']):
            util.stop_err('The file for the {} intervals ({}) does not exist.'.format(name, intrv_info['file']))

        instance = TargetCollection(intrv_info['file'])
        filtered_intervals = instance.filter_intervals(wanted_gene)
        saved_intervals = instance.save_intervals(filtered_intervals)
        if name == 'primary':
            primary_exons = saved_intervals
        intervals_list.append(saved_intervals)

    # Get union of the intervals and then merge
    combined_intervals = util.interval_union(*intervals_list)
    combined_intervals_merged = util.merge_intervals(combined_intervals, min_dist=min_dist)

    add_interval_labels(combined_intervals_merged, primary_exons)

    return combined_intervals_merged
