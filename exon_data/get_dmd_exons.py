from genepeeks.common import utilities as util
import pandas as pd

DMD_ensembl = util.Mongo.get_collection_data('gene', wanted_db='prod', query={'_id': 'DMD'}, find_one=True, single_field='ensembl')
DMD_exons = util.get_nested_value(DMD_ensembl, ('is_primary', 'transcripts', 'is_primary', 'exons'))
# DMD_exons_merged = util.merge_intervals(DMD_exons, min_dist=200)

exon_list = [['Ensembl_ID', 'Start', 'End']]
for exon in DMD_exons:
    exon_list.append([exon['id'], exon['start'], exon['end']])

df = pd.DataFrame(exon_list)
df.to_csv("exons.csv")
