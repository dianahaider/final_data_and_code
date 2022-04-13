#1.0 Running all trim lengths
for flen in {30..300..10} ##outer loop
do
	for rlen in {0..300..10} ##inner loop
	do
		mkdir -p "F"$flen"R"$rlen
		conda activate qiime2-2020.111 ##activate environment for trimming
		~/Documents/escuela/phd/plugin_paper/mock_code/18S/diana_18s/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-prok-scripts/P03-DADA2.sh $flen $rlen
		mv 03-DADA2d "F"$flen"R"$rlen
	done
done

#2.0 Pool all ASVs into one table
conda activate qiime2-2020.111

find ~/Documents/escuela/phd/plugin_paper/mock_code/16S/02-PROKs/alltrims -name 'table.qza' >all_tables.txt

python 02_a_consolidate_tables.py -i all_tables.txt -o merged_all_tables.tsv

echo 'Merged tables successfully'

#3.0 Assign taxonomy to each table
#3.1 find all the representative sequences from all trim lengths
find ~/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims -name 'representative_sequences.qza' >all_rep_seqs.txt
echo 'Found all representative sequences'

#3.2 run taxonomy for all of them
while read in; do
   qiime feature-classifier classify-sklearn \
    --i-reads "$in" \
    --i-classifier ~/Documents/escuela/phd/plugin_paper/mock_code/18S/silva-138-99-nb-classifier.qza \
    --output-dir "$in"-taxa --verbose; done < all_rep_seqs.txt
echo 'Generated all taxonomy'

#3.3 batch unzip de classification.qza and let result stay in parent folder
find . -name 'classification.qza' -exec sh -c 'unzip -d "${1%.*}" "$1"' _ {} \;

#rename all taxonomy.tsv to their trimlengths
python3 move_rename.py alltrims


mv -i F*R*.tsv all_taxonomies/ #puts all tsvs in new directory with correct names

python adapt_metadata(all_merged, manifest, metadata)

#add the
find ~/Documents/escuela/phd/plugin_paper/mock_code/16S/02-PROKs/alltrims/all_taxonomies -name 'F*R*.tsv' > all_taxonomies.txt

#run Bacaros at different levels

git clone https://github.com/alexmanuele/Bacaros_Beta.git

for i in {1..7..1}
do
	python run_beta.py --input ~/Documents/escuela/phd/plugin_paper/mock_code/16S/02-PROKs/alltrims/all_taxonomies.txt --metric t --l $i --output results-$i
done



conda deactivate



#
