conda activate qiime2-2020.111
print "QIIME2 activated"

find ~/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims -name 'table.qza' >all_tables.txt
print "Found all tables"

python consolidate_tables.py -i all_tables.txt -o merged_all_18s_tables.tsv

conda deactivate
