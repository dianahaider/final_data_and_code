#download alexmanuele/Bacaros_Beta from github
conda activate pluginanalysis

git clone https://github.com/alexmanuele/Bacaros_Beta.git

##batch unzip all classification qza files in their own filepath
find . -name 'classification.qza' -exec sh -c 'unzip -d "${1%.*}" "$1"' _ {} \;

python3 move_rename.py testfiles

#this makes all classifications renamed by their trim lengths and moved to a new directory

mv -i *.tsv taxonomies/ #puts all tsvs in new directory with correct names



#if need to batch change txt to tsv
#for f in *.txt; do mv -- "$f" "${f%.txt}.tsv"; done
