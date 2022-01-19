conda activate qiime2-2020.111

#find all the representative sequences from all trim lengths
find ~/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims -name 'representative_sequences.qza' >all_rep_seqs.txt
echo 'Found all representative sequences'

#run taxonomy for all of them
while read in; do
   qiime feature-classifier classify-sklearn \
    --i-reads "$in" \
    --i-classifier ~/Documents/escuela/phd/plugin_paper/mock_code/18S/silva-138-99-nb-classifier.qza \
    --output-dir "$in"-taxa --verbose; done < all_rep_seqs.txt
echo 'Generated all taxonomy'

find ~/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims -type f -name 'taxonomy.tsv' -exec cat {} + >all_18s_taxo.tsv

conda deactivate
