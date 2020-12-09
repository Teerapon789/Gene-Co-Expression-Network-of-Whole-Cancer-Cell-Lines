### title: "NetworkAnalysis"
### author: "Teerapon Sahwangarrom"

source("/Users/teeraponsahwangarrom/Desktop/RNASeq/FunctionFiles/networkFunctions.R")

data1 = loadAsList(
  "/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/Expression/030-OutlierSamplesRemoved/CancerRNASeqexpr-xdata-sampleAnnotation.OSR.RData")

expr = data1$expr.OSR
sampleAnnot = data1$sampleAnnotation.OSR
sampleAnnot.bin = sampleAnnot[, -1]

geneAnnot = data1$geneAnnotation.he

all.equal(geneAnnot$Name, colnames(expr))

rm(data1); gc()

#weights0 = loadTable(file = file.path("../Data/Expression/030-OutlierSamplesRemoved",
#            "CancerRNASeqxdataWeightFactorsForIOR.csv.gz"),
#          sep = ",", header = TRUE, transpose = TRUE, convertToMatrix = TRUE)

data2 = loadAsList(
  "/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/Expression/030-OutlierSamplesRemoved/CancerRNASeqxdataWeightsAndFactors.RData")

weights = data2$weightsForIOR[, match(colnames(expr), colnames(data2$weightsForIOR))]
#factors = data2$weightFactorsForIOR[, match(colnames(expr), colnames(data2$weightFactorsForIOR))]

rm (data2); gc();

## Scale Free Topology Analysis
enableWGCNAThreads(2)
powers = c(3:9)
sft = pickSoftThreshold(expr, dataIsExpr = TRUE, weights = weights, networkType = "signed hybrid", powerVector = powers, blockSize = 2000, verbose = 2)

R.sq = sft$fitIndices$SFT.R.sq
mean.k = sft$fitIndices$mean.k
median.k = sft$fitIndices$median.k

colors = "blue"
pdf(file = "Plots/scaleFreeTopologyAnalysis.pdf", wi = 12, he = 3)
#sizeGrWindow(12, 3)
par(mfrow = c(1,3))
par(mar = c(3.5, 3.5, 1.5, 1))
par(mgp = c(2.1, 0.8, 0))
plot(powers, R.sq, pch = 21, col = 1, bg = colors, cex = 1.7, main = "Scale-free topology fit index", xlab = "Soft-thresholding power", ylab = "Fit R squared")
addGrid(v = TRUE, linesPerTick = 2)

plot(powers, log10(mean.k), pch = 21, col = 1, bg = colors, cex = 1.7, main = "Mean connectivity", xlab = "Soft-thresholding power", ylab = "log10(Mean connectivity)")
addGrid(v = TRUE, linesPerTick = 2)

plot(powers, log10(median.k), pch = 21, col = 1, bg = colors, cex = 1.7, main = "Median connectivity", xlab = "Soft-thresholding power", ylab = "Median connectivity")
addGrid(v = TRUE, linesPerTick = 2)


## Network Construction
multiExpr = multiData(Cancer = expr)
multiWeights = multiData(Cancer = weights)
softPower = 4 

networkOptions = newNetworkOptions(correlationOptions = newCorrelationOptions(corType = "pearson", maxPOutliers = 0.05, corOptions = list(use = 'p', nThreads = 0)), power = softPower, networkType = "signed hybrid", TOMType = "signed", TOMDenom = "mean")

TOMDir.indiv= "/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/RData/TOMs"

if (!checkRDataFile("/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/RData/individualTOMinfo.RData", "individualTOMs"))
{
  dir.create(TOMDir.indiv, recursive = TRUE);
  gc();
  individualTOMs = individualTOMs(
    multiExpr, networkOptions = networkOptions,
    individualTOMFileNames = spaste(TOMDir.indiv, "/individualTOM-Set%s-Block%b.RData"),
    maxBlockSize = 45000,
    saveTOMs = TRUE,
    verbose = 5)
  
  save(individualTOMs, file = spaste("/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/RData/individualTOMinfo.RData"));
}

consensusTree = newConsensusTree( consensusOptions = newConsensusOptions(consensusQuantile = 0,
                                                                         calibration = "none"),
                                  inputs = "Cancer",
                                  analysisName = "Cancer");

if (!checkRDataFile("/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/RData/consensusTOMInfo.RData", "consensusTOMInfo"))
{
  system.time({ consensusTOMInfo = hierarchicalConsensusTOM(
    consensusTree = consensusTree,
    consensusTOMFilePattern = spaste("/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/RData/TOMs/consensusTOM-%a-Block%b.RData"),
    useDiskCache = FALSE, 
    individualTOMInfo = individualTOMs,
    keepIntermediateResults = FALSE,
    verbose = 10)});
  
  save(consensusTOMInfo, file = "/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/RData/consensusTOMInfo.RData")
}

gc();

mods.ns = hierarchicalConsensusModules(
  consensusTOMInfo = consensusTOMInfo,
  multiExpr = multiExpr,
  multiWeights = multiWeights, 
  checkMissingData = FALSE,
  detectCutHeight = 0.995,
  deepSplit = 2.5,
  minModuleSize = 30,
  useBranchEigennodeDissim = TRUE,
  mergeCutHeight = 0.20,
  minCoreKME = 0.4, minKMEtoStay = 0.2,
  iteratePruningAndMerging = TRUE,
  verbose = 5)

save(mods.ns, file = "/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/RData/mods.ns.RData")

labels = mods.ns$labels

geneSignificance = cor(expr, sampleAnnot.bin, weights.x = weights)


## Plot gene dendrogram, module colors and associations with cancers
moduleColors = labels2colors(labels)
colors.GS  = numbers2colors(geneSignificance, signed = TRUE, commonLim = TRUE)

plotColors = cbind(moduleColors, colors.GS)

plotTree = clipTree(mods.ns$dendrograms[[1]], p = 0.005)
rowLabels = spaste(colnames(geneSignificance), " (", sapply(sampleAnnot.bin, sum), ")")

pdf(file = "Plots/geneDendrogram.pdf", wi = 16, he = 9)
plotDendroAndColors(plotTree, plotColors, c("Modules", rowLabels), rowText = labels, textPositions = 1, cex.rowText = 1, rowWidths = c(1, 5, rep(1, ncol(geneSignificance))), rowTextAlignment = "center", dendroLabels = FALSE, main = spaste("Gene clustering, module colors and association with traits" ), marAll = c(0, 17, 2, 0), hang = 0.02, guideHang = 0.08, cex.colorLabels = 0.9, addGuide = TRUE)


## Output
MEs = moduleEigengenes(expr, labels, excludeGrey = TRUE)$eigengenes
rownames(MEs) = rownames(expr)
multiMEs = multiData(Cancer = MEs)

tab1 = hierarchicalConsensusKME(multiExpr, labels, multiWeights = multiWeights, multiEigengenes = multiMEs, setNames = names(multiExpr), consensusTree = consensusTree, useRankPvalue =FALSE, getOwnModuleZ = FALSE, getBestModuleZ = FALSE, additionalGeneInfo = geneAnnot[, -1, drop = FALSE],reportWeightType = "equal", includeWeightTypeInColnames = FALSE)

kmeTable = tab1[, !multiGrepl(c("meta.Z.", "kME.*Cancer"), names(tab1))]

names(kmeTable) = multiSub(c("consensus.", "[Cc]ons", "kME"), c("", "", "KME"), names(kmeTable))

dir.create("Results", recursive = TRUE)
write.csv(signifNumeric(kmeTable, 3), file = gzfile("/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/Results/networkAnalysisResults.csv.gz"), row.names  = FALSE, quote = TRUE)

write.csv(dataForTable(MEs, transpose = FALSE, IDcolName = "Sample"), file = "/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/Results/moduleEigengenes.csv", quote = FALSE, row.names = FALSE)

# Save samples

set.seed(123)

sam1 = sample(nrow(kmeTable), 100)
kmeTable.sample = kmeTable[sam1, ]

write.csv(signifNumeric(kmeTable.sample, 3), file = gzfile("/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/Results/networkAnalysisResults-sample.csv.gz"), row.names  = FALSE, quote = TRUE)

sam2 = sample(nrow(MEs), 50)
MEs.sample = MEs[sam2, ]

write.csv(dataForTable(MEs.sample, transpose = FALSE, IDcolName = "Sample"),file = "/Users/teeraponsahwangarrom/Desktop/GeneNetworkAnalysis/Results/moduleEigengenes-sample.csv",quote = FALSE, row.names = FALSE);


