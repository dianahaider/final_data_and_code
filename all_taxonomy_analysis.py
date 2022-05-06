#import libraries in python3 kernel
import pandas as pd
import seaborn as sns
import glob
import os
import boto
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import skbio
#import fastcluster #this thing is tricky for clustermap when not all relationship exists
from functools import reduce
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import GC
from collections import defaultdict
from collections import Counter
import statistics
import itertools as it
from scipy import stats
from matplotlib.ticker import FormatStrFormatter
%matplotlib inline

#paths to files needed as inputs for respective functions
all_merged_path = '/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/merged_all_18s_tables.tsv'
manifestfile_path = '/Users/Diana/Documents/escuela/phd/plugin_paper/final_data_and_code/input_data/MANIFEST.tsv'
metadatafile_path = '/Users/Diana/Documents/escuela/phd/plugin_paper/final_data_and_code/input_data/METADATA.tsv'
path_to_tsvs = '/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims/all_taxonomies/'
path_to_bacaro = '/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/TD_metric/Bacaros_Beta/all18sstagg_s_lvl7/all_18s_stagg_taxo.csv'
expected_file = 'expected_18s_staggered'
path_to_expected = '/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/IN-SILICO/in-silico-mocks/18s_stagg_even_expected_tax.tsv'
path_to_fastQC_forward = '/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/02-EUKs/run46/Forward/merged_f_fastqc/avg_base_q_forward.txt'
path_to_fastQC_reverse = '/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/02-EUKs/run46/Reverse/merged_r_fastqc/avg_base_q_reverse.txt'
path_to_best_taxo = '/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims/separated_taxonomies/Staggered/46F260R120.tsv'
path_to_compare = '/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims/separated_taxonomies/Staggered/46F20R0.tsv'
merged_fasta = '/Users/Diana/Documents/escuela/phd/plugin_paper/final_data_and_code/analyse_data/BLAST/merged.fasta'
blast_results = '/Users/Diana/Documents/escuela/phd/plugin_paper/final_data_and_code/analyse_data/BLAST/blast_results.txt'

paths_to_lvls = {
#    "Domain": ""
  "Phylum": "/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/TD_metric/Bacaros_Beta/all18slvl2/all_18s_stagg_taxo.csv",
  "Class": "/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/TD_metric/Bacaros_Beta/all18slvl3/all_18s_stagg_taxo.csv",
  "Order": "/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/TD_metric/Bacaros_Beta/all18slvl4/all_18s_stagg_taxo.csv",
  "Family": "/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/TD_metric/Bacaros_Beta/all18slvl5/all_18s_stagg_taxo.csv",
  "Genus": "/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/TD_metric/Bacaros_Beta/all_18s_stagg_lvl6/all_18s_stagg_taxo.csv",
  "Species": "/Users/Diana/Documents/escuela/phd/plugin_paper/mock_code/18S/TD_metric/Bacaros_Beta/all_18s_stagg_lvl7/all_18s_stagg_taxo.csv"
}



#2.0 Pool all ASVs into one table
conda activate qiime2-2020.111

find ~/Documents/escuela/phd/plugin_paper/mock_code/16S/02-PROKs/alltrims -name 'table.qza' >all_tables.txt

python 02_a_consolidate_tables.py -i all_tables.txt -o merged_all_tables.tsv

echo 'Merged tables successfully'


def divide_by_comm(all_merged, manifestfile, metadatafile, community, composition, runnumber):
    df = pd.read_csv(all_merged, sep='\t')
    tables = df[['sample_name', 'feature_id', 'feature_frequency']].copy()
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
    merged = pd.merge(merged,metadata, on='sample-id')
    merged = merged.replace({'V2': '16S'}, regex=True)

    print('Set up metadata ...')
    merged.to_csv('filtered_asvs.tsv', sep = '\t')
    print('Saved filtered_asvs.tsv')

#make df of features/composition+run+comm
    files = glob.glob('{0}*.tsv'.format('input_data/all_taxonomies_'+community+'/'))
    taxos = []
#    if not os.path.exists(path+composition):
#        os.mkdir(path+composition)
    for filename in files:
        tax = pd.read_csv(filename, sep='\t')
        tax['table_id'] = str(filename.split('/')[-1])
        tax["table_id"] = tax["table_id"].str.replace(".tsv", "")
        tax['Forward_trim'], tax['Reverse_trim'] = tax['table_id'].str.split('R', 1).str
        tax['Forward_trim'] = tax['Forward_trim'].map(lambda x: x.lstrip('F'))
        tax["Forward_trim"] = pd.to_numeric(tax["Forward_trim"])
        tax["Reverse_trim"] = pd.to_numeric(tax["Reverse_trim"])
        taxos.append(tax)
    print('Appended all taxonomies to taxos')
    taxos = pd.concat(taxos)
    taxos = taxos.rename(columns={"Feature ID": "feature_id"}, errors="raise")
    separated = taxos.merge(merged, how='left', on='feature_id')
    separated = separated.drop_duplicates()
    separated = separated[separated["community"] == community]
    separated = separated[separated["composition"] == composition]
    separated['run-number']= separated['run-number'].astype(str)
    separated = separated[separated["run-number"] == runnumber]
    separated['sum'] = separated.groupby(['table_id','sample-id'])['feature_frequency'].transform('sum')
    separated['ratio'] = separated['feature_frequency']/(separated['sum'])
    #make a dictionary with keys for id-ing the taxon belonging to this sub-community
    separated_dic = pd.Series(separated.Taxon.values,separated.feature_id.values).to_dict()

#generate folder of split taxonomies by runnumber and composition
    # Directory
    directory = composition+runnumber
    # Parent Directory path
    parent_dir = 'input_data/all_taxonomies_'+community
    # Path
    path = os.path.join(parent_dir, directory)
    # Create the directory
    # 'GeeksForGeeks' in
    # '/home / User / Documents'
    os.mkdir(path)
    for filename in files:
        taxonomy = pd.read_csv(filename, sep='\t')
        taxonomy = taxonomy.rename(columns={"Feature ID": "feature_id"}, errors="raise")
        newz = taxonomy.merge(merged, how='left', on='feature_id')
        #new = newz.drop(['sample-id'], axis=1)
        new = newz.drop_duplicates()
        new = new[new["community"] == community]
        new = new[new["composition"] == composition]
        new['run-number']= new['run-number'].astype(str)
        new = new[new["run-number"] == runnumber]
        new = new[new.feature_frequency != 0]
        new = new.rename(columns={"feature_id":"Feature ID"}, errors="raise")
        new = new[['Feature ID', 'Taxon', 'Confidence']].copy()
        new = new.drop_duplicates()
        new.to_csv(filename.split('all_taxonomies_'+community)[0]+'all_taxonomies_'+community+'/'+composition+runnumber+'/'+runnumber+filename.split('all_taxonomies_'+community+'/')[1], sep = '\t')
    return (files, separated, separated_dic, merged)


#get expstagg
## to import the expected taxonomies and transofrm to ratios
def expected_df(path):
    expected = pd.read_csv(path, sep='\t')
    expected = expected.rename(columns={'silva_taxonomy':'Taxon', 'sample-id': 'Replicate'})
    expected_even = expected[expected.mock_even_insilico != 0]
    expected_even = expected_even.drop(columns=['mock_staggered_insilico','taxonomy'])
    expected_even.reset_index(drop=True, inplace=True)
    expected_even['expected_ratio'] = expected_even['mock_even_insilico']/(expected_even['mock_even_insilico'].sum())
    expected_stagg = expected.drop(columns=['mock_even_insilico','taxonomy'])
    expected_stagg.reset_index(drop=True, inplace=True)
    expected_stagg['expected_ratio'] = expected_stagg['mock_staggered_insilico']/(expected_stagg['mock_staggered_insilico'].sum())
    expstaggplotting = expected_stagg.copy()
    expstaggplotting['family'] = [s.split('; ')[0] for s in expstaggplotting['group']]
    expstaggplotting['genus'] = [s.split('; ')[1] for s in expstaggplotting['group']]
    expstaggplotting['species'] = [s.split('; ')[2] for s in expstaggplotting['group']]

    return (expected_even, expected_stagg, expstaggplotting)



#run 1 sample t test
#takes trimcombination in format F40R40
def ttst(trimcombination, expdf):
    compr_obs_exp = expdf.merge(separated, how='outer', on='Taxon')
    compr_obs_exp = compr_obs_exp.rename(columns={'sample-id': 'Replicate'})
    compr_obs_exp=compr_obs_exp.fillna(0)
    compr_obs_exp["group"].replace({0: "False positive"}, inplace=True)
    df_with_groups = compr_obs_exp[compr_obs_exp.table_id == trimcombination]


    means = df_with_groups.groupby('group').mean()
    newsi=means[['expected_ratio', 'ratio']].copy()
    newsi.sort_values('expected_ratio', ascending = False)
    newsi['Difference'] = newsi['expected_ratio'] / newsi['ratio']
    newsi

    # generate a boxplot to see the data distribution by treatments. Using boxplot, we can easily detect the differences between different treatments
    ax = sns.boxplot(x='ratio', y='group', data=df_with_groups).set(
    xlabel='Relative abundance',
    ylabel='Group'
    )
    #ax.tick_params(axis='x', labelrotation=90)
    plt.show()

    groups = df_with_groups['group'].unique()

    results = []
    for group in groups:
        arr = df_with_groups[(df_with_groups['group'] == group)]['ratio'].values  #Filter the dataframe.
        results.append({'group': group,
                        'ratios': arr}) #Make a single "record" containing the table id, replicate, and ratio array.
        r_gr = pd.DataFrame.from_records(results)

    y = expdf['expected_ratio']
    for i in range(len(y[0:-1])):
        xi = r_gr.iloc[i]['ratios']
        yi = y[i]
        print(r_gr.iloc[i]['group'])
        print(stats.ttest_1samp(xi, yi))
        #but there are outliers
    return (r_gr, newsi)


#make r2 heatmap
def r2hm(expdf, outputname):
    grouped = compr_obs_exp.groupby(['Forward_trim', 'Reverse_trim'])
    test = (grouped.apply(lambda x: pd.Series(stats.linregress(x['expected_ratio'], x['ratio'])))
               .rename(columns={
                        0: 'slope',
                        1: 'intercept',
                        2: 'rvalue',
                        3: 'pvalue',
                        4: 'stderr'}))
    test['r2'] = test['rvalue']**2
    testcopy = test.reset_index()
    testcopy["Forward_trim"].replace({0: 280}, inplace=True)
    testcopy["Reverse_trim"].replace({0: 290}, inplace=True)
    tohm = testcopy.pivot("Forward_trim", "Reverse_trim", "r2")
    tohm.fillna(0)
    ax = sns.heatmap(tohm,cmap=("YlOrBr"))
    ax.invert_yaxis()
    plt.figure(figsize=(12,12))
    ax.set(xlabel='Reverse trim length', ylabel='Forward trim length')
    fig = ax.get_figure()
    fig.savefig(outputname+'.png', bbox_inches = "tight")

def get_fig_per_group(groupname, expectedratio, duplicate='mean'):
    neoceratium = separated[separated['Taxon'].str.contains(groupname)]
    neoceratium = neoceratium[neoceratium.feature_frequency !=0]
    neoceratium.rename(columns = {'sample-id':'sample_id'}, inplace=True)
    if duplicate!='mean':
        neoceratiumR1 = neoceratium[neoceratium.sample_id == 'R46-18S-'+duplicate]
    else:
        neoceratiumR1 = neoceratium.groupby(['Forward_trim','Reverse_trim'])[['ratio']].mean()
        neoceratiumR1 = neoceratiumR1.reset_index()
    neoceratiumR1["Forward_trim"] = pd.to_numeric(neoceratiumR1["Forward_trim"])
    neoceratiumR1["Reverse_trim"] = pd.to_numeric(neoceratiumR1["Reverse_trim"])
    neoceratiumR1["Forward_trim"].replace({0: 280}, inplace=True)
    neoceratiumR1["Reverse_trim"].replace({0: 290}, inplace=True)
    neoceratiumR1merged = neoceratiumR1.groupby(['Forward_trim','Reverse_trim'])[['ratio']].mean()
    neoceratiumR1merged = neoceratiumR1merged.reset_index()
    tohm = neoceratiumR1merged.pivot("Forward_trim", "Reverse_trim", "ratio")
    tohm.rename({280: 'full'}, axis=0, inplace=True)
    tohm.rename({290: 'full'}, axis=1, inplace=True)
    ax = sns.heatmap(tohm, cmap="Oranges_r")#, mask= (tohm < (expectedratio-(0.0005*expectedratio))) & (tohm > (expectedratio+(0.0005*expectedratio)))) #cmap=sns.color_palette("hls", 90)
    ax = sns.heatmap(tohm, mask=tohm <= expectedratio, cmap=sns.color_palette("GnBu", 5)) #square=True, annot=False, vmin=0, vmax=1, cbar=False, ax=ax)
    #ax = sns.heatmap(tohm, mask=tohm >= expectedratio, cmap=sns.color_palette("Oranges_r", 5)) #square=True, annot=False, vmin=0, vmax=1, cbar=False, ax=ax)
    ax.invert_yaxis()
    plt.figure(figsize=(12,12))
    ax.set(xlabel='Reverse trim length', ylabel='Forward trim length')
    fig = ax.get_figure()
    fig.savefig('ratio'+groupname+'.png', bbox_inches = "tight")
    return neoceratiumR1merged

def blast_analysis(blastoutput):
    #upload the dataframe and give headers
    df= pd.read_csv(blastoutput,
                sep='\t',
               names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    pd.to_numeric(df['evalue'])
    blastoutput = df.copy()
    #rename query sequence id column (observed feature id)
    blastoutput = blastoutput.rename(columns={'qseqid': 'feature_id'})
    #persid = balstoutput.groupby(['sseqid', 'feature_id']).mean().reset_index()
    #reassign the names
    blastoutput = blastoutput.merge(id_keys,on='sseqid')
    blastoutput = blastoutput[blastoutput.qstart == 1]
    blastoutput = blastoutput[blastoutput.length > 250]
    merged = blastoutput.merge(separated, how='outer', on='feature_id')
    significant_only = merged[merged['evalue'] <1e-50 ]
    groupedmerged = merged.groupby(['short_id','feature_id'])[['evalue']].agg('min') #was max
    groupedmerged = groupedmerged.reset_index()
    groupedmerged = groupedmerged.dropna()
    tohm = groupedmerged.pivot("short_id", "feature_id", "evalue")
    tohm=tohm.fillna(0)

    #plot
    kws = dict(cbar_kws=dict(orientation='horizontal'))
    ax = sns.clustermap(tohm, metric="euclidean", norm=LogNorm(), standard_scale=1, method="ward", cmap=sns.color_palette("mako", 200),
                   mask=(tohm==0), **kws)
    ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xmajorticklabels())#, fontsize = 5)

    threshold = 1e-200
    x_labels_ticks = ax.ax_heatmap.get_xticklabels()
    total_asvs_above_threshold = 0
    for i, xtickdata in enumerate(x_labels_ticks):
        asv = xtickdata._text
        if tohm[asv].max() <= threshold:
            total_asvs_above_threshold = total_asvs_above_threshold + 1
        else:
            xtickdata._text = ''

    # reset the tick labels with the modified list
    ax.ax_heatmap.set_xticklabels(x_labels_ticks)
    hm = ax.ax_heatmap
    hm.set_xlabel("Observed ASVs")
    hm.set_ylabel("Expected ASVs")
    x0, _y0, _w, _h = ax.cbar_pos
    ax.ax_cbar.set_position([x0, 0.9, ax.ax_row_dendrogram.get_position().width, 0.02])
    ax.ax_cbar.tick_params(axis='x')
    plt.show()
    ax.fig.savefig('blast_evalues.png', format='png', dpi=600, bbox_inches="tight")

    return merged, groupedmerged


def make_fasta(pathmergedfasta, outputfilename, R='separated', F='separated'):
    if R!='separated':
        rallfs = separated[separated.Reverse_trim == R]
        separated_dic = pd.Series(rallfs.Taxon.values,rallfs.feature_id.values).to_dict()
    else:
        separated_dic = pd.Series(separated.Taxon.values, separated.feature_id.values).to_dict()
    if F!='separated':
        fallrs = separated[separated.Forward_trim == F]
        separated_dic = pd.Series(fallrs.Taxon.values,fallrs.feature_id.values).to_dict()
    else:
        separated_dic = pd.Series(separated.Taxon.values, separated.feature_id.values).to_dict()

    fa = SeqIO.parse(pathmergedfasta,
                 "fasta")
    seqs_i_want = [] #we'll put the good sequences here
    for record in fa: #a SeqRecord has the accession as record.id, usually.
        if record.id in separated_dic.keys(): #This is how you check if the accession is in the values of the dict
            seqs_i_want.append(record)
    #Now we can write the list of records to a fasta file. This will take care of the formatting etc
    with open(outputfilename+'.fasta', "w") as f:
        SeqIO.write(seqs_i_want, f, "fasta")
    return print('Saved selected sequences as '+outputfilename+'.fasta')

def hm(path, truthfilename):
    bacaros_dm = pd.read_csv(path)
    bacaros_dm = bacaros_dm.set_index('Unnamed: 0')
   # bacaros_dm = 1  - bacaros_dm
    #bacaros_dm is a distance matrix of table X table
    #my_pcoa = skbio.stats.ordination.pcoa(bacaros_dm.values)
    #plt.scatter(my_pcoa.samples['PC1'],  my_pcoa.samples['PC2'])
    against_exp = bacaros_dm[[truthfilename]].copy()
    against_exp = against_exp.reset_index().rename(columns={against_exp.index.name:'sample_name'})
    against_exp.drop(against_exp.index[against_exp['sample_name'] == truthfilename], inplace=True)
    against_exp['Forward_trim'] = [s.split('R')[0] for s in against_exp['sample_name']]
    against_exp['Forward_trim'] = [s.split('46F')[1] for s in against_exp['Forward_trim']]
    against_exp['Reverse_trim'] = [s.split('R')[1] for s in against_exp['sample_name']]
    against_exp["Forward_trim"] = pd.to_numeric(against_exp["Forward_trim"])
    against_exp["Reverse_trim"] = pd.to_numeric(against_exp["Reverse_trim"])
    against_exp["Forward_trim"].replace({0: 280}, inplace=True)
    against_exp["Reverse_trim"].replace({0: 290}, inplace=True)
    tohm = against_exp.pivot("Forward_trim", "Reverse_trim", "expected_18s_staggered")
    tohm.rename({280: 'full'}, axis=0, inplace=True)
    tohm.rename({290: 'full'}, axis=1, inplace=True)
    ax = sns.heatmap(tohm, cmap=sns.color_palette("YlOrBr", 200), vmin=0, vmax=1) #cmap=sns.color_palette("hls", 90)

    ax.invert_yaxis()
    plt.figure(figsize=(12,12))
    ax.set(xlabel='Reverse trim length', ylabel='Forward trim length')
    fig = ax.get_figure()
    fig.savefig('bacaros_deltaT.png', bbox_inches = "tight")

    return (tohm, bacaros_dm, against_exp)


def get_thresholds():
    thresholds = np.arange(0, 1, 0.05).tolist()
    df = pd.DataFrame(columns=['Rank', *thresholds]) #* unpack the list and adds to list, **is for unpacking dictionaries
    #eg. d3 = {**d1, **d2} or l3 = [*l1, *l2]
    for key in dic.keys():
        tohm, bacaros_dm, agexp = hm(dic[key], expected_file) #input the value of each key from the dictionary
        row ={'Rank': key} #make an empty dictionary
        for thresh in thresholds:
            pp = ((tohm[tohm > thresh ].count().sum())/812)*100
            row[thresh] = pp
        df = df.append(row, ignore_index=True)
    newdf = df.set_index(['Rank'])
    newdf['Rank'] = newdf.index
    melted = newdf.melt('Rank', var_name='Thresholds', value_name='TD')
    melted[["Thresholds", "TD"]]=melted[["Thresholds", "TD"]].apply(pd.to_numeric)
    sns.set_context('paper')
    fig = plt.figure(figsize=(5,5))
    p = sns.lineplot(data=melted, x="Thresholds", y="TD", hue='Rank')
    p.set_xlabel("\u0394T") #, fontsize = 10)
    p.set_ylabel("Percentage of tables",)# fontsize = 10)
    #p.invert_xaxis()
    ax = p.get_figure()
    ax.savefig('thresholds.png', format='png', dpi=600, bbox_inches="tight")
    return df

    import sys
    if __name__ == "__main__":
        # parse arguments
