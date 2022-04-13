#import libraries in python3 kernel
import pandas as pd
import glob
import os
import numpy as np


def split_tx_by_sample(path, composition, runnumber, foldername, community):
    files = glob.glob('{0}*.tsv'.format(path))
    for filename in files:
        tax = pd.read_csv(filename, sep='\t')
        tax = tax.rename(columns={"Feature ID": "feature_id"}, errors="raise")
        new = tax.merge(merged, how='left', on='feature_id')
        new = new.drop(['Confidence', 'sample-id'], axis=1)
        new = new.drop_duplicates()
        new = new[new["community"] == community]
        new = new[new["composition"] == composition]
        new['run-number']= new['run-number'].astype(str)
        new = new[new["run-number"] == runnumber]
        new.to_csv(filename.split(foldername)[0]+foldername+'/'+composition+'/'+runnumber+filename.split(foldername+'/')[1], sep = '\t')
    return new
