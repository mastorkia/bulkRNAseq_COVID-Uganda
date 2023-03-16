# bulk-RNA-seq scripts
These scripts were used for the analysis of 100 bulk-RNA-seq samples studied in the "Multidimensional host profiling of COVID-19 in Uganda reveals prognostic immunometabolic signatures that persist during HIV co-infection and diverge by variant-driven epidemic phase" paper.

### Pipeline steps
* Adaptor removal and quality trimming using Trimmomatic.
* Reads mapping and gene count by STAR.
* Filtering, normalization and differential expression steps based on DESeq2 package.