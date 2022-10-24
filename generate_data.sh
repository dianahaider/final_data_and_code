#0.0 Download data and clone git repos with pipeline
conda env create -f environment.yml

git clone https://github.com/jcmcnch/eASV-pipeline-for-515Y-926R.git
#change the qiime version from the pipeline to your current qiime2 version
find . -type f | xargs perl -pi -e 's/qiime2-2019.4/qiime2-2020.11/g'

chmod a+x eASV-pipeline-for-515Y-926R/./ #makes files executable
source eASV-pipeline-for-515Y-926R/DADA2-pipeline/00-trimming-sorting-scripts/00-run-cutadapt.sh
conda env create -f eASV-pipeline-for-515Y-926R/env/bbmap.yaml #removed some packages that failed
source eASV-pipeline-for-515Y-926R/DADA2-pipeline/00-trimming-sorting-scripts/01-sort-16S-18S-bbsplit.sh #make sure you change the paths for the in-silico, and for the bbsplit db, this step is long

#1.0 Running all trim lengths
#FOR PROKARYOTES
cd 02-PROKs
source ~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-prok-scripts/P00-create-manifest.sh
source ~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-prok-scripts/P01-import.sh

for flen in {30..300..10} ##outer loop
do
	for rlen in {0..300..10} ##inner loop
	do
		mkdir -p all_trims/"F"$flen"R"$rlen
		conda activate qiime2-2020.111 ##activate environment for trimming
		~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-prok-scripts/P03-DADA2.sh $flen $rlen
		mv 03-DADA2d all_trims/"F"$flen"R"$rlen
	done
done
echo 'Generated all 16S trim combinations'

#FOR EUKARYOTES
cd ~/MOCK_ANALYSIS/02-EUKs
source ~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E00-create-manifest-viz.sh
source ~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E01-import.sh

for flen in {260..300..20} ##outer loop
do
	for rlen in {260..300..20} ##inner loop
	do
		mkdir -p all_trims/"F"$flen"R"$rlen
		conda activate bbmap-env ##activate environment for trimming
		~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E03-bbduk-cut-reads.sh $flen $rlen
		mv 03-size-selected all_trims/"F"$flen"R"$rlen
		cd all_trims/"F"$flen"R"$rlen
		~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E04-fuse-EUKs-withoutNs.sh
		~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E05-create-manifest-concat.sh
		conda activate qiime2-2020.111
		~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E06-import-concat.sh
		~/MOCK_ANALYSIS/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E08-DADA2.sh
		cd ..
		cd ..
	done
done
echo 'Generated all 18S trim combinations'
cd ~/MOCK_ANALYSIS

#2.0 Assign taxonomy to each table
#2.1 find all the representative sequences from all trim lengths
find . -name 'representative_sequences.qza' >all_rep_seqs.txt
echo 'Found both 16S and 18S representative sequences'

#2.2 Download qiime2 pre-fitted on SILVA 99% naive bayes classifier, or download personalized classifier
# If on Mac OS and need to install wget, uncomment the next line
# brew install wget
wget https://data.qiime2.org/2022.2/common/silva-138-99-nb-classifier.qza

#2.3 Run taxonomic assignment for all of them
#this should also be ran in your version of q2
while read in; do
   qiime feature-classifier classify-sklearn \
    --i-reads "$in" \
    --i-classifier silva-138-99-nb-classifier.qza \
    --output-dir "$in"-taxa --verbose; done < all_rep_seqs.txt
rm all_rep_seqs.txt
echo 'Generated all taxonomy'


#2.4 batch unzip the classification.qza and leave contents in parent folder
find . -name 'classification.qza' -exec sh -c 'unzip -d "${1%.*}" "$1"' _ {} \;

#3.0 Run in-silico sequences

#import in silico separately to get the right taxonomy? dig into the taxonomy files; have to run this
cd in-silico-mock
ls -d */ | sed 's#/##' >all_comms.txt #only folders without the slash

while read in; do
	cd $in
	qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path MANIFEST.tsv --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred64V2
	qiime vsearch join-pairs --i-demultiplexed-seqs paired-end-demux.qza --o-joined-sequences joined_sequences.qza;
	qiime vsearch dereplicate-sequences --i-sequences joined_sequences.qza --o-dereplicated-sequences dereplicated-sequences.qza --o-dereplicated-table dereplicated-table.qza
	qiime feature-classifier classify-sklearn --i-reads dereplicated-sequences.qza --i-classifier ~/MOCK_ANALYSIS/silva-138-99-nb-classifier.qza --output-dir taxonomy
	unzip taxonomy/classification.qza
	cd ..; done < all_comms.txt

rm all_comms.txt

#4.0 Batch unzip all DADA2 denoising_stats.qza and leave contents in parent folder
find . -name 'denoising_stats.qza' -exec sh -c 'unzip -d "${1%.*}" "$1"' _ {} \;

#5.0 Extract average nucleotide score per read through QIIME2 
#mkdir for each community (in ~/Diana/data/)
#add manifest file in each folder
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path pe-64-manifest --output-path reads.qza --input-format PairedEndFastqManifestPhred33
qiime demux summarize --i-data reads.qza --o-visualization demux.qzv
#open or batch unzip them and reach for the forward-seven-number-summaries.tsv in data folder
