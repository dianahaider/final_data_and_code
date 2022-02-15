#tables are merged to metadata based on sample_id
#first get all tables from 02_get_table, outputs merged_all_18s_tables.tsv
#then run these
#need to take as input the merged all tables

#import libraries in python3 kernel
import pandas as pd
import glob
import boto
import numpy as np
import skbio
import statistics
from pathlib import Path
import sys
import os

#path is the path to the folder that contains the taxonomy .tsv files
#composition and runnumber are the metadata variable you want
#foldername is the name (not path) of folder that contains the taxonomies
def load_df(path, composition, runnumber, adaptedmetadata):
    files = glob.glob('{0}*.tsv'.format(path))
#    if not os.path.exists(path+composition):
#        os.mkdir(path+composition)
    Path(path+composition).mkdir(parents=True, exist=_ok)
    for filename in files:
        merged = pd.read_csv(adaptedmetadata, sep='\t')
        tax = pd.read_csv(filename, sep='\t')
        tax = tax.rename(columns={"Feature ID": "feature_id"}, errors="raise")
        new = tax.merge(merged, how='left', on='feature_id')
        new = new.drop(['sample-id'], axis=1)
        new = new.drop_duplicates()
        new = new[new["community"] == '18S']
        new = new[new["composition"] == composition]
        new['run-number']= new['run-number'].astype(str)
        new = new[new["run-number"] == runnumber]
        new = new[['feature_id', 'Taxon', 'Confidence']].copy()
        new.to_csv(path+composition+'/'+runnumber+filename.split('/')[-1], sep = '\t')

if __name__ == '__main__':
    load_df(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
