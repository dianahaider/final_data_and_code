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

def adapt_metadata(all_merged, manifestfile, metadatafile, outdir):
    df = pd.read_csv(all_merged, sep='\t')
    tables = df[['sample_name', 'feature_id']].copy()
    tables.rename(columns={'sample_name':'file'}, inplace=True)
    manifest = pd.read_csv(manifestfile, sep='\t')
    manifest['file'] = [s.split('SPOT_USC_2/')[1] for s in manifest['absolute-filepath']]
    manifest['file'] = [s.split('.R')[0] for s in manifest['file']]
    manifest = manifest.drop(columns = ['absolute-filepath', 'direction'])
    manifest.drop_duplicates()
    merged = pd.merge(tables,manifest, on='file')
    merged = merged.drop(columns = ['file'])
    merged = merged.drop_duplicates()
    print('Set up manifest ...')
    metadata = pd.read_csv(metadatafile, sep='\t')
    metadata.rename(columns={'sample-id':'sample_name'}, inplace=True)
    merged = pd.merge(tables,metadata, on='sample_name')
    merged = merged.rename(columns={"run-number": "run_number"}, errors="raise")
    merged = pd.merge(merged,metadata, on='sample-id')
    merged.to_csv(outdir+'filtering_asvs.tsv', sep = '\t')



if __name__ == '__main__':
    adapt_metadata(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv)
