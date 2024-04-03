# bulk-RNA-seq scripts
These scripts were used for the analysis of 100 bulk-RNA-seq samples studied in the "COVID-19 immune signatures in Uganda persist in HIV co-infection and diverge by pandemic phase" paper [PMID:38368384](https://pubmed.ncbi.nlm.nih.gov/38368384/).

### Pipeline steps
* Adaptor removal and quality trimming using Trimmomatic.
* Reads mapping and gene count by STAR.
* Filtering, normalization and differential expression steps based on DESeq2 package.
