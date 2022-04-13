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

id_keys = pd.DataFrame(np.array([['in.silico.stag_1', 'Acineta_flava_KR-10010701_Acinetidae_X_sp._strain5', 'Ciliophora; Phyllopharyngea; Suctoria'],
                       ['in.silico.stag_81', 'Amoebophrya_sp._Amoebophrya_sp.', 'Dinophyta; Syndiniales; Dino-Group-II_a'],
                       ['in.silico.stag_721', 'Amoebophrya_sp._Dino-Group-II-Clade-10-and-11_X_sp.', 'Dinophyta; Syndiniales; Dino-Group-II_b'],
                       ['in.silico.stag_1041', 'Ceratium_longipes_ccmp1770_Ceratium_tenue', 'Dinophyta; Dinophyceae; Neoceratium'],
                       ['in.silico.stag_1043', 'Chrysochromulina_simplex_partial_Chrysochromulina_X_sp.', 'Haptophyta; Prymnesiophyceae; Prymnesiales'],
                       ['in.silico.stag_1723', 'Cryothecomonas_sp._APCC_Cryothecomonas_sp.', 'Cercozoa; Filosa-Thecofilosea; Cryomonadida'],
                       ['in.silico.stag_2043', 'Guinardia_delicatula_Guinardia_delicatula', 'Ochrophyta; Bacillariophyta; Coscinodiscophyceae_a'],
                       ['in.silico.stag_2523', 'Gymnodinium_sp._Dinophyceae', 'Dinophyta; Dinophyceae; Gyrodinium'],
                       ['in.silico.stag_2803', 'Gymnodinium_sp._Gymnodinium_dorsalisulcum','Dinophyta; Dinophyceae; Gymnodinium'],
                       ['in.silico.stag_3083', 'Larcopyle_butschlii_Larcopyle_butschlii','Radiolaria; RAD-B; Larcopyle'],
                       ['in.silico.stag_3283', 'Leptocylindrus_convexus_SZN-B768_Radial-centric-basal-Coscinodiscophyceae_X_sp.','Ochrophyta; Bacillariophyta; Coscinodiscophyceae_b'],
                       ['in.silico.stag_3337', 'Lingulodinium_polyedrum_Lingulodinium_polyedrum','Dinophyta; Dinophyceae; Lingulodinium'],
                       ['in.silico.stag_3737', 'Paracalanus_parvus_Paracalanus_parvus','Metazoa; Arthropoda; Crustacea'],
                       ['in.silico.stag_3741', 'Strombidium_cf._basimorphum_Strombidium_basimorphum','Ciliophora; Spirotrichea; Oligotrichia'],
                       ['in.silico.stag_3781', 'Woloszynskia_cincta_Woloszynskia_cincta','Dinophyta; Dinophyceae; Dinophyceae_X'],
                       ['in.silico.stag_3801', 'Karlodinium_micrum_clone_Dino-Group-I-Clade-5_X_sp.','Dinophyta; Syndiniales; Dino-Group-I'],
                      ]),
             columns=['sseqid', 'taxonomic_id', 'short_id'])


#input alltables merged
def adapt_metadata(all_merged, manifestfile, metadatafile):
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
    print('Set up metadata ...')
    merged.to_csv('filtering_asvs.tsv', sep = '\t')
    print('Saved filtering_asvs.tsv')
    return merged


#get separated
def get_abundances(path, composition, runnumber):
    files = glob.glob('{0}*.tsv'.format(path))
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
    separated = separated[separated["community"] == '18S']
    separated = separated[separated["composition"] == composition]
    separated['run-number']= separated['run-number'].astype(str)
    separated = separated[separated["run-number"] == runnumber]
    separated['sum'] = separated.groupby(['table_id','sample-id'])['feature_frequency'].transform('sum')
    separated['ratio'] = separated['feature_frequency']/(separated['sum'])
    #make a dictionary with keys for id-ing the taxon belonging to this sub-community
    separated_dic = pd.Series(separated.Taxon.values,separated.feature_id.values).to_dict()

    return (separated, separated_dic)


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
