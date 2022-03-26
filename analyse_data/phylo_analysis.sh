#find and unzip all fasta files
find ~/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims -name 'representative_sequences.qza' -exec unzip {} \;

#concat all sequences from all trims (18S only)
cat ~/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims/*/data/dna-sequences.fasta > merged.fasta

#moved the file to final_data_and_code/analysis/phylogenetics
#move to jup notebook now and run function select_fasta, which uses dictionary made from get_abundances
#move output to final_data_and_code/analysis/phylogenetics folder

#remove duplicates
awk '/^>/{f=!d[$1];d[$1]=1}f' outputidi.fasta > dedup_18s.fa

#make phylotree
FastTree -gtr -nt < 04_alignedmafft.fasta > tree_18stg.nwk



#go to jup notebook
