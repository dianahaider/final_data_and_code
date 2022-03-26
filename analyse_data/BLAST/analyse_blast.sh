#make the merged sequences file of observed
find ~/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims -name 'representative_sequences.qza' -exec unzip {} \;
cat ~/Documents/escuela/phd/plugin_paper/mock_code/18S/all_trims/*/data/dna-sequences.fasta > merged.fasta

#install blast from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#Create the subject database
#Since it's the expected sequence, we don't need to inspect the quality plots, we can directly fuse euks into one read using bbmap-env fuse
#!/bin/sh
bash E04-fuse-EUKs-withoutNs.sh
#Use the concatenated reads to make a db
#Since each sequence is present many times to represent different proportions in the community, extract only the unique sequences
#this 18STAGG_18s_mock_staggered_insilico_Rconcatenated.fastq contains 4000 sequences of 16 clones in different proportions
cat 18STAGG_18s_mock_staggered_insilico_Rconcatenated.fastq | paste - - - - | LC_ALL=C sort -t$'\t' -k2,2 -u | tr "\t" "\n" > unique18SSTAG.fq
cat unique18SSTAG.fq | paste - - - - |cut -f 1, 2| sed 's/@/>/'g | tr -s "/t" "/n" > unique18SSTAG.fa
#can double check if the correct #sequences by 'wc -l out.fq' and divide by 4
#Make it into a blast db

makeblastdb -in unique18SSTAG.fa -out unique18SSTAG_db.fa -dbtype 'nucl' -hash_index -parse_seqids
#run blast against all diana_unique_sequences
#Make query database
blastn -db ~/Documents/escuela/phd/plugin_paper/mock_code/IN-SILICO/in-silico-mocks/04-concatenated/unique18SSTAG_db.fa -query ~/Documents/escuela/phd/plugin_paper/mock_code/IN-SILICO/in-silico-mocks/04-concatenated/merged.fasta -outfmt 6 -out ~/Documents/escuela/phd/plugin_paper/mock_code/IN-SILICO/in-silico-mocks/04-concatenated/blast_results.txt -num_threads 4
#continue analysis on python
#install seqkit to parse through fastq headers
