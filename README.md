# GeneNetworkAnalysis

This analysis I had developed the R code for preprocessing RNA-seq data and network analysis, which I used the knowledge from the fundamental of WGCNA (Langfelder and Horvath 2008).

The script pre-processes the raw counts by filtering out low-expressed features, running an approximate variance-stabilizing transformation, and calculating observation whose purpose is to downweigh outliers. Also, the script creates rudimentary sample traits by inferring the type of cancer from sample names and creating binary indicator variables for cancers that have at least five samples.

The network analysis includes scale free topology analysis, network construction, plotting gene dendrogram, module colours and association with cancers, and functional enrichment analysis of the module genes (GO terms).
