#tables are merged to metadata based on sample_id
#first get all tables from 02_get_table, outputs merged_all_18s_tables.tsv
#then run these
#need to take as input the merged all tables

def load_df(filenames, composition, runnumber, foldername):
    for filename in filenames:
        tax = pd.read_csv(filename, sep='\t')
        tax = tax.rename(columns={"Feature ID": "feature_id"}, errors="raise")
        new = tax.merge(merged, how='left', on='feature_id')
        new = new.drop(['Confidence', 'sample-id'], axis=1)
        new = new.drop_duplicates()
        new = new[new["community"] == '18S']
        new = new[new["composition"] == composition]
        new['run_number']= new['run_number'].astype(str)
        new = new[new["run_number"] == runnumber]
        new.to_csv(filename.split(foldername)[0]+foldername+'/'+composition+'/'+runnumber+filename.split(foldername+'/')[1], sep = '\t')
