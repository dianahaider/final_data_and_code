# 515F/926R mock communities analysis script

### Structure of the starting repository

```Working directory  
└── Cloned repo  
    ├── 00-raw  
    │   ├── raw_reads.fastq.gz  
    ├── in-silico-mocks  
    │   ├── 02-EUKs/02-PROKs  
    │          └── Even/Staggered  
    │                 ├── reads.fastq  
    ├── generate_data.sh
    ├── run_analysis.ipynb
    ├── environment.yml
    ├── MANIFEST.tsv
    └── METADATA.tsv
 ```
    
### Datasets    
| Dataset       | Community     | Number of ASVs     |
| ------------- | ------------- | -------- |
| 16S           | Even          | NewYork  |
| 18S           | Test2         | Toronto  |

### Generate the data
Run



```Working directory  
└── Cloned repo  
    ├── 00-raw  
    │   ├── raw_reads.fastq.gz  
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
    │                 ├── reads.fastq  
    │                 └── taxonomy  
    │                      └── classification.qza     
    ├── generate_data.sh
    ├── run_analysis.ipynb
    ├── environment.yml
    ├── MANIFEST.tsv
    └── METADATA.tsv
```
    

