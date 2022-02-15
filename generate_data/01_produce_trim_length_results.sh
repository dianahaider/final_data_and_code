for flen in {0..300..10} ##outer loop
do
	for rlen in {0..300..10} ##inner loop
	do
		mkdir -p "F"$flen"R"$rlen
		conda activate bbmap-env ##activate environment for trimming
		~/mocks-to-share-with-Diana/diana_18s/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E03-bbduk-cut-reads.sh $flen $rlen
		mv 03-size-selected "F"$flen"R"$rlen
		cd "F"$flen"R"$rlen
		~/mocks-to-share-with-Diana/diana_18s/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E04-fuse-EUKs-withoutNs.sh
		~/mocks-to-share-with-Diana/diana_18s/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E05-create-manifest-concat.sh
		conda activate qiime2-2020.111
		~/mocks-to-share-with-Diana/diana_18s/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E06-import-concat.sh
		~/mocks-to-share-with-Diana/diana_18s/eASV-pipeline-for-515Y-926R/DADA2-pipeline/01-euk-scripts/E08-DADA2.sh
		cd ..
	done
done
