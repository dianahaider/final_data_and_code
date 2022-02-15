#To obtain all the tables, the pipeline from https://www.protocols.io/view/default-fuhrman-lab-pipeline-for-exact-amplicon-se-vi9e4h6 is run without changing parameters other than the trim length, therefore the manifest and metadata files need to be adapted to the manually curated metadata file to split ASVs by communities and run numbers

#import libraries in python3 kernel
import pandas as pd
import glob
import boto
import numpy as np
import skbio
import statistics
import sys

#takes the path to each file

df = pd.read_csv('/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/merged_all_18s_tables.tsv', sep='\t')
#Inspect the columns to keep only relevant ones
tables = df[['sample_name', 'feature_id', 'feature_frequency']].copy()
#modify the sample names to matching ones with metadata
tables.rename(columns={'sample_name':'file'}, inplace=True)
manifest['file'] = [s.split('SPOT_USC_2/')[1] for s in manifest['absolute-filepath']]
manifest['file'] = [s.split('.R')[0] for s in manifest['file']]
#drop unused columns
manifest = manifest.drop(columns = ['absolute-filepath', 'direction'])
merged = pd.merge(tables,manifest, on='file')
merged = merged.drop(columns = ['file'])
merged = merged.drop_duplicates()
#the sample name is what will allow me to merge with the metadata
metadata = pd.read_csv('/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/16S/Runs_46_47_55_18years_SPOT_USC_2/METADATA_F08.tsv', sep='\t')
#Merge tables on sample name
merged = pd.merge(merged,metadata, on='sample-id')

def get_abundances(path, composition, runnumber):
    files = glob.glob('{0}*.tsv'.format(path))
    taxos = []
#    if not os.path.exists(path+composition):
#        os.mkdir(path+composition)
    for filename in files:
        tax = pd.read_csv(filename, sep='\t')
        tax['table_id'] = str(filename.split('/')[-1])
        tax['Forward_trim'] = str(s.split('R')[0] for s in tax['table_id'])
        tax['Forward_trim'] = str(s.split('46F')[1] for s in tax['Forward_trim'] )
        tax['Reverse_trim'] = str(s.split('R')[1] for s in tax['table_id'])
        taxos.append(tax)
    taxos = pd.concat(taxos)
    taxos = taxos.rename(columns={"Feature ID": "feature_id"}, errors="raise")
    new = taxos.merge(merged, how='left', on='feature_id')
    return new
