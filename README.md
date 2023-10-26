# 515F/926R mock communities analysis script

We determined the stability of ASV clusters in terms of taxonomy, abundance and phylogeny generated with DADA2 through a variety of trimming tresholds.

### Structure of the starting repository

```Working directory  
└── Cloned repo  
    ├── 00-raw  
    │   └── raw_reads.fastq.gz  
    ├── in-silico-mocks  
    │   └── 02-EUKs/02-PROKs  
    │          └── Even/Staggered  
                      ├── clone_names.tsv
    │                 └── reads.fastq  
    ├── generate_data.sh
    ├── run_analysis.ipynb
    ├── environment.yml
    ├── MANIFEST.tsv
    └── METADATA.tsv
 ```
    
### Datasets
| Dataset       | Community     | Number of clones   | Link | 
| ------------- | ------------- | ------------------ |------|
| 16S           | Even          | 11                 |https://www.ebi.ac.uk/ena/browser/view/PRJEB48162
| 16S           | Staggered     | 27                 |https://www.ebi.ac.uk/ena/browser/view/PRJEB48162
| 18S           | Even          | 10                 |https://www.ebi.ac.uk/ena/browser/view/PRJEB35673
| 18S           | Staggered     | 16                 |https://www.ebi.ac.uk/ena/browser/view/PRJEB35673

### Generate the data
Run ```generate_data.sh``` from the cloned repo, then open the ```run_analysis.ipynb``` jupyter notebook and run all.

### Final repository

```Working directory  
└── Cloned repo  
    ├── 00-raw  
    │   └── raw_reads.fastq.gz  
    ├── 01-trimmed  
    │   └── trimmed_reads.fastq  
    ├── 02-EUKs/02-PROKs  
    │   ├── 00-fastq  
    │   │   └── cdhit_euk_reads.fastq  
    │   ├── intermediate_files  
    │   │   ├── all_taxonomies  
    │   │   │   ├── Even/Staggered
    │   │   │   │     └── Taxonomies.tsv  
    │   │   └── all_trims  
    │   │       └── Forward and reverse trim combinations  
    │   │             └── 03-DADA2d  
    │   │                 ├── denoising_statistics.qza  
    │   │                 ├── representative_sequences.qza  
    │   │                 ├── table.qza  
    │   │                 └── taxonomy  
    │   │                      └── classification.qza  
    │   └── logs  
    ├── in-silico-mocks  
    │   └── 02-EUKs/02-PROKs  
    │          └── Even/Staggered  
    |                 ├── clone_names.tsv
    │                 ├── reads.fastq  
    │                 └── taxonomy  
    │                      └── classification.qza     
    ├── generate_data.sh
    ├── run_analysis.ipynb
    ├── environment.yml
    ├── MANIFEST.tsv
    └── METADATA.tsv
```
