### title: "GeneNetworkAnalysis"
###author: "Teerapon Sahwangarrom"

source("/Users/teeraponsahwangarrom/Desktop/RNASeq/FunctionFiles/networkFunctions.R")
source("/Users/teeraponsahwangarrom/Desktop/RNASeq/FunctionFiles/preprocessing.R")
source("/Users/teeraponsahwangarrom/Desktop/RNASeq/FunctionFiles/outlierRemovalFunctions.R")

library(DESeq2)


data = read.table("/Users/teeraponsahwangarrom/Desktop/CCLE_RNAseq_081117.reads.txt.bz2", sep="\t", header = TRUE)

GeneAnnot0 = data[, 1:2]
counts = t(as.matrix(data[, -c(1:2)]))
colnames(counts) = GeneAnnot0$Name

type = sub("[^_]*_", "", rownames(counts))
table(type)
typeDF = data.frame(type = type)
type.bin = binarizeCategoricalColumns(typeDF, val1 = 0, val2 = 1, includePrefix = FALSE, minCount = 5, includePairwise = FALSE, includeLevelVsAll = TRUE, dropFirstLevelVsAll = FALSE, dropUninformative = TRUE, nameSep = "")

sampleAnnot = cbind(typeDF, type.bin)


ds = DESeqDataSetFromMatrix.PL(t(counts), sampleAnnot, design = ~1)

ds = estimateSizeFactors(ds)

counts.norm = t(counts(ds, normalized = TRUE))
zz = heColumns(counts, 8, 0.1)

x = log2(counts.norm[, zz] + 1)
var = colVars(x)
means = colMeans(x)

av = approximateVST(counts.norm[, zz])


prepr = preprocessGeneral( xdata = counts, addToAnalysisManager = FALSE, organism = "human", normalizationPipeline = "DESeq2", idsAreSymbols = FALSE, analysisName = "CancerRNASeq", pheno = sampleAnnot, sampleAnnotation = sampleAnnot, plotPheno = sampleAnnot, geneAnnotation = GeneAnnot0, minValue = 10, minProportion = 0.1, useApproximateVST = TRUE, vst.design = "~1", bw.groupBy = "type", bw.minSamplesPerGroup = 1, outlierRemovalZ = 8, ior.replace = FALSE, ior.remove = FALSE, saveXData.entrez = FALSE, saveDir.base = "/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis", stp.mainBase = "Cancer RNA Seq", fileNamePrefix = "CancerRNASeq", verbose = 4)

