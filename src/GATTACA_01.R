# Header Info ----------------------------------------------------------------------------------------------------------
#
# GATTACA v0.1-alpha - General Algorithm for The Transcription Analysis by one-Channel Arrays
#
# a FeAR R-script - 13-Nov-2020
#
# Pipeline for one-Color (HD) Microarrays
# Data are supposed to be already background-subtracted, log2-transformed, and interarray-normalized
#





# * Package Loading ----------------------------------------------------------------------------------------------------

library(preprocessCore) # For Quantile-Quantile Normalization
library(rafalib)        # For Bland Altman Plot (aka MA Plot) 
library(PCAtools)       # Principal Component Analysis
library(genefilter)     # Expression Gene Filtering
library(limma)          # Empirical Bayes Method for Differential Expression 
library(RankProd)       # Rank Product Method for Differential Expression
library(VennDiagram)    # Generate Venn Diagrams
library(openxlsx)       # Read, Write, and Edit .xlsx (Excel) Files





# Prepare Annotations --------------------------------------------------------------------------------------------------

system.root = "D:\\UniUPO Drive\\" # SilverLife @ Home
system.root = "D:\\Drive UniUPO\\" # SkyLake2   @ DBIOS

# Affymetrix Human Genome U133 Set (A) annotation data (chip hgu133a)
library(hgu133a.db) 
#hgu133a() # List of the available annotations
annot = data.frame(Accession   = sapply(contents(hgu133aACCNUM), paste, collapse = ", "),
                   GeneSymbol  = sapply(contents(hgu133aSYMBOL), paste, collapse = ", "),
                   Description = sapply(contents(hgu133aGENENAME), paste, collapse = ", "))

# Affymetrix Human Genome U133 Set (A) Local Annotation File
setwd(paste(system.root, "Coding\\R scripts\\Annotations\\HG-U133A-na36-annot-csv", sep = ""))
annot = read.xlsx("HG-U133A.na36.annot.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
#as.matrix(colnames(annot)) # List of the available annotations
annot = annot[,c("Representative_Public_ID", "Gene_Symbol", "Gene_Title")]

# Affymetrix Human Genome U133 Set (B) annotation data (chip hgu133b)
library(hgu133b.db)
#hgu133b() # List of the available annotations
annot = data.frame(Accession   = sapply(contents(hgu133bACCNUM), paste, collapse = ", "),
                   GeneSymbol  = sapply(contents(hgu133bSYMBOL), paste, collapse = ", "),
                   Description = sapply(contents(hgu133bGENENAME), paste, collapse = ", "))

# Affymetrix Human Genome U133 Set (B) Local Annotation File
setwd(paste(system.root, "Coding\\R scripts\\Annotations\\HG-U133B-na36-annot-csv", sep = ""))
annot = read.xlsx("HG-U133B.na36.annot.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
#as.matrix(colnames(annot)) # List of the available annotations
annot = annot[,c("Representative_Public_ID", "Gene_Symbol", "Gene_Title")]

# Agilent SurePrint G3 Mouse GE 8x60K Local Annotation File
setwd(paste(system.root, "Coding\\R scripts\\Annotations\\AllAnnotations - Agilent Update 2018", sep = ""))
annot = read.xlsx("028005_D_AA_20181026.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
#as.matrix(colnames(annot)) # List of the available annotations
annot = annot[,c("GeneSymbol", "GeneName", "Description")]





# * Variable Definition ------------------------------------------------------------------------------------------------
# User-Defined Experiment-Specific Variables

myFolder = "C:\\Users\\aleph\\Desktop"
myFolder = "D:\\Dropbox\\temp for today"
myFile = "1 - log_Intensity_matrix_Organized_all.txt"

rowOffset = 1    # Row offset (rows to skip, including the header)
colWithNames = 1 # Column containing (unique) gene identifiers
colOffset = 1    # Column offset (columns to skip, including row names)

groups = c("WT","Ab","AbFK","TG","TGFK") # Experimental Design - Group Names
design = c(5,5,5,5,4) # Experimental Design - Ordered Group Size (compact Design Mode)

# ...or build a custom Experimental Design vector (Full Design Mode)
# Use ALL natural numbers from 1 to m=length(groups) - According to 'groups' ordering
design = c(3,3,1,1,1,1,1,4,4,4,5,5,2,2,2,2,2,3,3,3,4,4,5,5)

myColors = c("cornflowerblue","firebrick3","olivedrab3","darkgoldenrod1","purple","magenta3")
if (length(myColors) < length(groups)) {
  cat("\n")
  stop("Too few colors in \'myColors\' vector!\n\n")
}

# Filter: log2-expression conventional threshold (to be checked from case to case)
thr0 = 6 # Agilent
thr0 = 4 # Affymetrix

# Fold Change Threshold. Usually, |log2FC|>0.5 OR |log2FC|>1
thrFC = 0.5

# If-statement control flags
saveOut = TRUE
logConversion = FALSE
secondNorm = FALSE
save.PNG.plot = TRUE
save.PDF.plot = FALSE
append.annot = TRUE





# * Function Definition ------------------------------------------------------------------------------------------------
# User-Defined Functions

# Save Graphical Output to 'GATTACA Figures' Sub-directory
printPlots = function(figureName, PNG.bool = save.PNG.plot, PDF.bool = save.PDF.plot)
{
  figSubFolder = "GATTACA Figures"
  fullName = file.path(figSubFolder, figureName, fsep = .Platform$file.sep) # OS-independent path separator
  
  if (!file.exists(figSubFolder) & (PNG.bool | PDF.bool)) {
    dir.create(figSubFolder)
    cat("\nNew folder '", figSubFolder, "' has been created in the PWD...\n\n", sep = "")
  }
  if (PNG.bool) { # invisible(capture.output()) to suppress automatic output to console
    invisible(capture.output(
      dev.print(device = png, filename = paste(fullName, ".png", sep = ""), width = 820, height = 600)))
  }
  if (PDF.bool) {
    invisible(capture.output(
      dev.print(device = pdf, paste(fullName, ".pdf", sep = ""))))
  }
}





# * Data Loading -------------------------------------------------------------------------------------------------------
# Gene Expression Matrix - log2-Intensity-Values

setwd(myFolder)

# Read raw data
dataset = read.table(myFile, header = FALSE, sep = "\t", dec = ".")
d = dim(dataset)
cat("\nRaw dataset dimensions:", d, "\n\n", sep = " ")
dataset[1:10,1:5]

# Extract column headings (sample identifiers)
header  = read.table(myFile, nrows = 1, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
dataset = read.table(myFile, skip = rowOffset, header = FALSE, sep = "\t", dec = ".")
colnames(dataset) = unlist(header)
d = dim(dataset)
cat("\nHeaded dataset dimensions:", d, "\n\n", sep = " ")
dataset[1:10,1:5]

# Extract row names (gene identifiers)
rownames(dataset) = dataset[,colWithNames]
dataset = dataset[,(colOffset+1):d[2]]
header  = header[,(colOffset+1):d[2]]
d = dim(dataset)
cat("\nFinal dataset dimensions:", d, "\n\n", sep = " ")
dataset[1:10,1:5]





# * Rename Samples  ----------------------------------------------------------------------------------------------------
# Tidy Sample Names According to the Experimental Design

sampleName = vector() # To declare an empty vector
m = length(groups)
badMsg = "Bad Experimental \'design\' vector!\n\n"

# Design check
if (length(design) != m && length(design) != d[2]) {
  cat("\n")
  stop(badMsg)
}
if (length(design) == m) {
  design = rep(c(1:m), design) # Convert from Compact to Full Design mode
  cat("\nCompact \'design\' mode approved\n\n")
} else if (length(design) == d[2]) {
  if (max(design) > m || min(design) < 1 || sum((design %% 1) != 0)) { # %% is R modulus operator (for integer checking)
    cat("\n")
    stop(badMsg)
  } else {
    cat("\nFull \'design\' mode approved\n\n")
  }
}

# Create a new vector containing tidy group names
for (i in 1:m) {
  index = which(design == i)
  grpName = paste(groups[i], "_", c(1:length(index)), sep = "") # Concatenate Strings with no separator
  sampleName[index] = grpName
}

# Correspondences Table
corrTable = cbind(unlist(header), sampleName) # Cast to matrix
colnames(corrTable) = c("Original_ID", "R_Name")
rownames(corrTable) = c(1:d[2])
corrTable
if (saveOut) {
  write.table(corrTable, "Corresp Table.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  cat("\n'Corresp Table.txt' has been saved in", myFolder, "\n\n")
}

colnames(dataset) = sampleName
dataset[1:10,1:5]





# Log Transformation ---------------------------------------------------------------------------------------------------
# Distribution Inspection to Spot Non-Logarithmic Data

boxplot(dataset, las = 2) # las = style of axis labels
plotDensities(dataset, legend = FALSE) # From limma package

# log2-transform iff intensities are in linear scale
if (logConversion) {
  
  # Are there any negative values?
  negVal = length(which(dataset < 0))
  if (negVal > 0) {
    cat("\n")
    stop(negVal, " negative expression values detected!\n\n")
  }
  
  # Are there any zero values?
  zeroVal = length(which(dataset == 0))
  cat("\n", zeroVal, "null values found\n\n")
  
  # Value imputation
  if (zeroVal > 0) {
    # Non-zero minimum value
    imputation = min(dataset[which(dataset > 0, arr.ind = TRUE)])
    #imputation = NaN # Alternatively
    
    dataset[which(dataset == 0, arr.ind = TRUE)] = imputation
    cat("\n", zeroVal, "null values have been imputated with", imputation, "\n\n")
  }
  
  length(which(dataset == 0))
  min(dataset)
  
  dataset = log2(dataset)
}





# Normalization --------------------------------------------------------------------------------------------------------
# After-RMA 2nd Quantile Normalization

if (secondNorm) {
  normData = normalize.quantiles(as.matrix(dataset)) # From preprocessCore package
  normData = as.data.frame(normData)
  rownames(normData) = rownames(dataset)
  colnames(normData) = colnames(dataset)
  
  # Inspect...
  normData[1:10,1:5]
  boxplot(normData, las = 2) # las = style of axis labels
  plotDensities(normData, legend = FALSE) # From limma package
  
  # ...if ok, then reassign
  dataset = normData
}





# MA-Plot & Box-Plot ---------------------------------------------------------------------------------------------------
# Normalization Final Check with Figure Production

boxplot(dataset, las = 2, col = myColors[design],
        main = "Expression values per sample", ylab = "log2 (intesity)") # las = style of axis labels
printPlots("1 - Boxplots")
plotDensities(dataset, legend = FALSE, main = "Expression values per sample") # From limma package
printPlots("2 - Density")

# MA-Plot for bias detection
plotMD(dataset, column = 1) # From limma package: array 'column' vs the average of all the arrays other than that
for (i in 1:(m-1)) { # All the possible combinations of two groups
  for (j in (i+1):m) {
    
    arr1 = rowMeans(dataset[,which(design == i)])
    arr2 = rowMeans(dataset[,which(design == j)])
    
    # From rafalib package
    maplot(arr1, arr2, 
           xlab = "A (Average log-expression)",
           ylab = "M (Expression log-ratio)",
           n = 5e4,
           curve.add = TRUE, curve.col = myColors[2], curve.lwd = 1.5, curve.n = 1e4,
           pch = 20, cex = 0.1)
    
    # As an alternative, from limma package (without trend curve)
    #plotMD(cbind(arr2, arr1),
    #       xlab = "A (Average log-expression)",
    #       ylab = "M (Expression log-ratio)")
    
    title(main = paste(groups[j], "vs", groups[i]))
    abline(h = 0, col = myColors[1], lty = 2) # lty = line type
    abline(h = c(1,-1), col = myColors[1])
    abline(v = thr0, col = myColors[1]) # Platform-specific log2-expression threshold
    printPlots(paste("3 - MA-Plot ", groups[j], "_vs_", groups[i], sep = ""))
  }
}





# Clustering -----------------------------------------------------------------------------------------------------------
# Sample-wise Hierarchical Clustering for Batch-Effect Detection

sampleDist = dist(t(dataset), method = "euclidean") # distance matrix: compute the distances between the ROWS of a matrix
hc = hclust(sampleDist, method = "ward.D")
plot(hc) # Display Dendrogram
printPlots("4 - Dendrogram")
# Desired number of clusters
kNum = 6

clust = cutree(hc, k = kNum) # Cut tree into kNum clusters
clustList = list() # Create an empty list
for (i in 1:kNum) {
  clustList[[i]] = cbind(clust[which(clust == i)]) # list of matrices
}
clustList

rect.hclust(hc, k = kNum, border = myColors[2]) # Red borders around the kNum clusters 
printPlots("4 - Dendrogram and Clusters")





# PCA ------------------------------------------------------------------------------------------------------------------
# Performed on Samples for Batch-Effect Detection

myMetadata = data.frame(row.names = colnames(dataset)) # Strictly enforced that rownames(metadata) == colnames(dataset)

# Add column 'Group' to metadata dataframe ($ allows accessing/creating dataframe columns by name)
myMetadata$Group = rep(NA, d[2])
for (i in 1:m) {
  myMetadata$Group[which(design == i)] = groups[i]
}

# Do the PCA (centering the data before performing PCA, by default)
pcaOut = pca(dataset, metadata = myMetadata)

# Plot results
screeplot(pcaOut)
printPlots("5 - Scree Plot")

colMAtch = myColors[1:m] # Vector for color matching
names(colMAtch) = groups
biplot(pcaOut, colby = "Group", colkey = colMAtch)
printPlots("6 - PCA")
pairsplot(pcaOut, colby = "Group", colkey = colMAtch)
printPlots("6 - PCA Pairs")

# Possibly remove some 'batched' sample, by sample name, e.g.: toBeRemoved = c("TG_1","WT_2","TGFK_1","Ab_5")
toBeRemoved = c()
if (length(toBeRemoved) > 0) {
  for (i in 1:length(toBeRemoved)) {
    rem.Index = which(colnames(dataset) == toBeRemoved[i])
    dataset = dataset[,-rem.Index]
    design = design[-rem.Index]
    sampleName = sampleName[-rem.Index]
  }
  cat("\n", length(toBeRemoved), " groups have been removed\n", sep = "")
  d = dim(dataset)
  cat("\nSub-dataset dimensions:", d, "\n\n", sep = " ")
  dataset[1:2,]
}





# SD vs Mean Plot ------------------------------------------------------------------------------------------------------
# Poisson Hypothesis Check

# Un-log intensity values
unlogged.dataset = 2^dataset
unlogged.dataset[1:10,1:5]

# Store values using matrices
meansArr = matrix(nrow = d[1], ncol = m+1)
SDsArr = matrix(nrow = d[1], ncol = m+1)
corrArr = vector()

# Statistics for each group...
for (i in 1:m) {
  meansArr[,i] = rowMeans(unlogged.dataset[,which(design == i)], na.rm = TRUE) # Within-group mean
  SDsArr[,i] = apply(unlogged.dataset[,which(design == i)], 1, sd, na.rm = TRUE) # Within-group SD
  corrArr[i] = cor(meansArr[,i], SDsArr[,i]) # Mean-SD Correlation
}
# ...and for the whole experiment
meansArr[,m+1] = rowMeans(unlogged.dataset, na.rm = TRUE)
SDsArr[,m+1] = apply(unlogged.dataset, 1, sd, na.rm = TRUE)
corrArr[m+1] = cor(meansArr[,m+1], SDsArr[,m+1])

# Scatter plot
par(mfrow = c(1, m+1)) # Optional, for sub-plotting
X.max = max(meansArr)
Y.max = max(SDsArr)
for (i in 1:(m+1)) {
  
  plot(meansArr[,i], SDsArr[,i],
       xlab = "Mean", ylab = "SD",
       xlim = c(0, X.max), ylim = c(0, Y.max),
       pch = 20, cex = 0.1)
  
  if (i <= m) {
    title(main = groups[i])
  }
  else {
    title(main = "Global")
  }
  mtext(side = 3, paste("Corr =", toString(round(corrArr[i], digits = 5))))
}
printPlots("7 - SD_vs_Mean Plot") # Save just the 'Global' one
par(mfrow = c(1, 1)) # To reset sub-plotting





# Filtering ------------------------------------------------------------------------------------------------------------
# Eliminate Low-Intensity Probes

kFac = 0.75 # Minimum gene presence per group - Default=0.80
grSize = vector()
for (i in 1:m) {
  grSize[i] = sum(design == i)
}
kSize = ceiling(grSize * kFac) # 'ceiling' to be more stringent than 'round'

# Filtering Table
filTable = cbind(grSize, grSize*kFac, kSize, round(100*(kSize/grSize), 1)) # Cast to matrix
colnames(filTable) = c("Group_Size", "Min_Presence", "Rounded", "Actual %")
rownames(filTable) = groups
cat("\nMinimum gene presence per group = ", kFac*100, "% of the samples\n\n", sep = "")
filTable
if (saveOut) {
  write.table(filTable, "Filtering Table.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  cat("\n'Filtering Table.txt' has been saved in", myFolder, "\n\n")
}

# Above thr0 value in at least kSize samples out of grSize (kFac*100 % - Default=80%)...
presenceIndx = matrix(nrow = d[1], ncol = m) # Store values using matrices
for (i in 1:m) {
  f1 = kOverA(kSize, thr0)
  ffun1 = filterfun(f1)
  presenceIndx[,i] = genefilter(dataset[,which(design == i)], ffun1)
}
retained = cbind(colSums(presenceIndx)) # Cast to matrix
colnames(retained) = c("Retained_Genes")
rownames(retained) = groups
retained

#...in at least 1 group (row-wise logical OR)
unionSet = (rowSums(presenceIndx) > 0)

filtSize = sum(unionSet)
cat("\nFiltering Report:\n",
    "-----------------\n",
    "  ", d[1] - filtSize, " genes have been filtered out (",
    round(((d[1] - filtSize)/d[1])*1e2, digits = 2), "%)\n",
    "  ", filtSize, " retained genes out of ", d[1], " (",
    round((filtSize/d[1])*1e2, digits = 2), "%)\n\n", sep = "")

dataset = dataset[which(unionSet),]





# * Contrast Definition ------------------------------------------------------------------------------------------------

filtSize = dim(dataset)[1] # Reassignment to get 'filtSize' defined even in the case of no filtering

counter = 1
myContr = vector()
for (i in 1:(m-1)) { # All the possible combinations of two groups
  for (j in (i+1):m) {
    myContr[counter] = paste(groups[j], "-", groups[i], sep = "") # Character vector specifying contrasts
    counter = counter + 1
  }
}

myContr = cbind(myContr) # Cast to matrix
rownames(myContr) = c(1:length(myContr))
myContr
# Modify here to retain (and reorder) just the contrasts of interest
myContr = myContr[c(3,6,1,2)]
myContr = cbind(myContr)
rownames(myContr) = c(1:length(myContr))
myContr





# * DE by Limma --------------------------------------------------------------------------------------------------------
# Differential Expression Assessment and Figure Production

limmaDesign = matrix(data = 0, nrow = d[2], ncol = m)
for (i in 1:m) {
  limmaDesign[which(design == i),i] = 1
}
colnames(limmaDesign) = groups
limmaDesign

fit = lmFit(dataset, limmaDesign)

contrast.matrix = makeContrasts(
  contrasts = myContr,
  levels = limmaDesign)
contrast.matrix

fit2 = contrasts.fit(fit, contrast.matrix)
efit2 = eBayes(fit2)

# Print Results (Top-Ten genes) for all the contrasts of interest
for (i in 1:length(myContr)) {
  # 'print' because automatic printing is turned off in loops (and functions)
  cat("\nDEG Top-List for contrast: ", myContr[i], "\n", sep = "")
  print(topTable(efit2, coef = i, adjust.method = "BH", sort.by = "B")) # just on-Screen
}

# Compute full DEG Tables
DEGs.limma = list() # Create an empty list
for (i in 1:length(myContr)) {
  DEGs.limma[[i]] = topTable(efit2, coef = i, number = filtSize,
                             adjust.method = "BH", sort.by = "B") # list of Data Frames
}

# Save full DEG Tables
if (saveOut) {
  for (i in 1:length(myContr)) {
    degTabName = paste("Limma - DEG Table ", myContr[i], ".txt", sep = "")
    write.table(DEGs.limma[[i]], degTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    cat("\n'", degTabName, "' has been saved in ", myFolder, "\n\n", sep = "")
  }
}

# Summary of DEGs (you can change the Log2-Fold-Change Threshold lfc...)
results.limma = decideTests(efit2, adjust.method = "BH", p.value = 0.05, lfc = thrFC)
summary(results.limma)
vennDiagram(results.limma)
printPlots("8 - Limma Venn")

# Show Hyperparameters
d0 = efit2$df.prior           # prior degrees of freedom
dg = mean(fit$df.residual)    # original degrees of freedom
hyp = cbind(c(efit2$s2.prior, # prior variance
      mean(fit$sigma^2),      # mean sample residual variance
      mean(efit2$s2.post),    # mean posterior residual variance
      d0, dg, d0/(d0+dg)))    # Shrinkage degree
rownames(hyp) = c("Prior Var",
                  "Sample Residual Mean Var",
                  "Posterior Residual Mean Var",
                  "Prior df",
                  "Original df",
                  "Shrinkage degree")
colnames(hyp) = "Hyperparameters"
hyp

# Save significant DEG list in Excel format with annotations
all.limma.sigs = list() # Create an empty list
for (i in 1:length(myContr)) {
  all.limma.sigs[[i]] = topTable(efit2, coef = i, number = filtSize,
                                 adjust.method = "BH", sort.by = "B",
                                 p.value = 0.05, lfc = thrFC) # list of Data Frames
}
for (i in 1:length(myContr)) {
  if (saveOut & dim(all.limma.sigs[[i]])[1] > 0) {
    if (append.annot) {
      # To merge two data frames horizontally (right outer join)
      joined = merge(annot, all.limma.sigs[[i]], by.x = "row.names", by.y = "row.names", all.y = TRUE)
      rownames(joined) = joined[,1]
      all.limma.sigs[[i]] = joined[,-1]
    }
    write.xlsx(all.limma.sigs[[i]], paste("Significant Genes by limma - ", myContr[i], ".xlsx", sep = ""),
               colNames = TRUE, rowNames = TRUE, sheetName = myContr[i],
               keepNA = TRUE, firstRow = TRUE) # Freezes the first row!
  }
}





# Limma Plot -----------------------------------------------------------------------------------------------------------
# MA-Plots with significant DEGs and Volcano plots

# Find Axis Limits
max.M.value = 0
for (i in 1:length(myContr)) {
  temp = max(abs(DEGs.limma[[i]][,1]))
  if (temp > max.M.value) {
    max.M.value = temp
  }
}
max.A.value = 0
for (i in 1:length(myContr)) {
  temp = max(DEGs.limma[[i]][,2])
  if (temp > max.A.value) {
    max.A.value = temp
  }
}
min.A.value = Inf
for (i in 1:length(myContr)) {
  temp = min(DEGs.limma[[i]][,2])
  if (temp < min.A.value) {
    min.A.value = temp
  }
}
min.P.value = 1
for (i in 1:length(myContr)) {
  temp = min(DEGs.limma[[i]][,4])
  if (temp < min.P.value) {
    min.P.value = temp
  }
}

# MA-Plot with DEGs
for (i in 1:length(myContr)) {
  # Mark in red/blue all the up-/down- regulated genes (+1/-1 in 'results.limma' matrix)
  plotMD(efit2, coef = i, status = results.limma[,i], values = c(1,-1), hl.col = myColors[c(2,1)],
         xlim = c(min.A.value, max.A.value), ylim = c(-max.M.value, max.M.value),
         xlab = "A (Average log-expression)",
         ylab = "M (log2-Fold-Change)")
  abline(h = 0, col = myColors[1], lty = 2) # lty = line type
  abline(h = c(thrFC,-thrFC), col = myColors[1])
  abline(v = thr0, col = myColors[1]) # Platform-specific log2-expression threshold
  printPlots(paste("9 - MA-Plot with Limma DEGs ", myContr[i], sep = ""))
}

# Volcano Plots - from limma package
for (i in 1:length(myContr)) {
  
  tot.DEG = sum(DEGs.limma[[i]][,5] < 0.05) # Total number of significant DEGs (without any FC cutoff)
  high.DEG = min(c(5,tot.DEG)) # To highlight no more than 5 genes per plot
  
  volcanoplot(efit2, coef = i, style = "p-value",
              highlight = high.DEG, names = rownames(dataset), hl.col = myColors[1],
              xlim = c(-max.M.value, max.M.value), ylim = c(0, -log10(min.P.value)),
              xlab = "log2-Fold-Change")
  title(main = myContr[i])
  
  # Significance Threshold
  # Find that p-value corresponding to BH-adj.p-value ~ 0.05 (or Bonferroni point when tot.DEG = 0)
  thrP = (0.05/filtSize)*(tot.DEG + 1)
  
  # Check the threshold p-value
  cat("\n", myContr[i], " - Threshold Report:\n",
      "--------------------------------------\n",
      "  p-value threshold  =  ", thrP, "\n",
      "  -log10(p-value)    =  ", -log10(thrP), "\n",
      "  Gene Ranking       =  ", tot.DEG, ":", tot.DEG + 1, "\n\n", sep = "")
  
  # Check in the DEG list
  # It should be the alpha-crossing-point for BH-adj.p-values column...
  # ...and a thrP-containing interval for un-adjusted p-values column.
  if (tot.DEG > 0) {
    print(DEGs.limma[[i]][tot.DEG:(tot.DEG + 1),])
    cat("\n")
  } else {
    print(DEGs.limma[[i]][1,])
    cat("\n")
  }
  
  if (FALSE) {
    # Alternative approach - Suitable also for correction methods other than BH
    indexP = which.min(abs(DEGs.limma[[i]][,5] - 0.05)) # Find the gene with adj.p.value closest to 0.05
    if (tot.DEG > 0) {
      thrP = DEGs.limma[[i]][indexP,4] # This is the p-value of the last DEG (or first non-DEG)
    } else {
      thrP = (0.05/filtSize) # Bonferroni point
    }
    
    # Check the threshold p-value
    cat("\nThreshold Report:\n",
        "-----------------\n",
        "  p-value threshold  =  ", thrP, "\n",
        "  -log10(p-value)    =  ", -log10(thrP), "\n",
        "  Gene Ranking       =  ", indexP, "\n\n", sep = "")
    
    # Check in the DEG list
    DEGs.limma[[i]][indexP:(tot.DEG + 1),]
  }
  
  abline(v = c(thrFC,-thrFC), col = myColors[2])
  abline(h = -log10(thrP), col = myColors[2])
  abline(h = -log10(0.05), col = myColors[2], lty = 2)
  printPlots(paste("10 - Volcano with Limma DEGs ", myContr[i], sep = ""))
}





# * DE by RankProduct --------------------------------------------------------------------------------------------------
# Differential Expression Assessment and Figure Production

results.RP = matrix(data = 0, nrow = filtSize, ncol = length(myContr))
rownames(results.RP) = rownames(dataset)
colnames(results.RP) = myContr
for (i in 1:length(myContr)) {
  
  contr.groups = strsplit(myContr[i], split = "-", fixed = TRUE)
  
  case.id = which(groups == contr.groups[[1]][1])
  case.index = which(design == case.id)
  ctrl.id = which(groups == contr.groups[[1]][2])
  ctrl.index = which(design == ctrl.id)
  
  sub.dataset = dataset[,c(ctrl.index,case.index)]
  cl = c(rep(0,length(ctrl.index)), rep(1,length(case.index)))
  
  # invisible(capture.output()) is to suppress automatic output to console
  # WARNING: therein <- (instead of =) is mandatory for assignment!
  invisible(capture.output(RP.out <- RankProducts(sub.dataset, cl, gene.names = rownames(dataset),
                                                 logged = TRUE, na.rm = FALSE, plot = FALSE, rand = 123)))
  invisible(capture.output(plotRP(RP.out, cutoff = 0.05)))
  
  # Compute full DEG Tables
  invisible(capture.output(DEGs.RP <- topGene(RP.out, logged = TRUE, logbase = 2, num.gene = filtSize)))
  for (j in 1:2) {
    DEGs.RP[[j]][,3] = log2(1/DEGs.RP[[j]][,3]) # Invert FC to get Case vs Ctrl and take the log2 values
    colnames(DEGs.RP[[j]])[3] = "Log2FC" # Correct column name
  }
  
  # Print Results (Top-Ten genes) for all the contrasts of interest
  cat("\nDEG Top-List for contrast: ", myContr[i], "\n", sep = "")
  tops = rbind(DEGs.RP$Table1[1:10,], DEGs.RP$Table2[1:10,])
  tops = tops[order(tops[,4])[1:10],] # Sort by PFP (~FDR) and take just the 'absolute' Top-Ten
  print(tops) # just on-Screen
  
  # Save full DEG Tables
  if (saveOut) {
    upDegTabName = paste("RP_Up - DEG Table ", myContr[i], ".txt", sep = "")
    dwnDegTabName = paste("RP_Down - DEG Table ", myContr[i], ".txt", sep = "")
    write.table(DEGs.RP[[1]], upDegTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    write.table(DEGs.RP[[2]], dwnDegTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    cat("\n'", upDegTabName, "' and '", dwnDegTabName, "' \nhave been saved in ", myFolder, "\n\n", sep = "")
  }
  
  # Fetch indexes of significant DEGs
  ups.Index = which(DEGs.RP$Table1[,4] < 0.05 & DEGs.RP$Table1[,3] > thrFC)
  dwn.Index = which(DEGs.RP$Table2[,4] < 0.05 & DEGs.RP$Table2[,3] < -thrFC)
  
  # Here 'results.RP' is created in analogy to 'results.limma': (-1,0,+1) for (Down, NotSig, Up)-regulated genes
  # Accessing matrix rows by row name: matrix[rowname, col_i]
  results.RP[rownames(DEGs.RP$Table1)[ups.Index],i] = 1
  results.RP[rownames(DEGs.RP$Table2)[dwn.Index],i] = -1
  
  # Save significant DEG list in Excel format with annotations
  if (saveOut & (length(ups.Index) > 1 | length(dwn.Index) > 1)) {
    
    all.RP.sigs = rbind(DEGs.RP$Table1[ups.Index,], DEGs.RP$Table2[dwn.Index,]) # To join two data frames vertically
    
    if (append.annot) {
      # To merge two data frames horizontally (right outer join)
      joined = merge(annot, all.RP.sigs, by.x = "row.names", by.y = "row.names", all.y = TRUE)
      rownames(joined) = joined[,1]
      all.RP.sigs = joined[,-1]
    }
    
    write.xlsx(all.RP.sigs, paste("Significant Genes by RP - ", myContr[i], ".xlsx", sep = ""),
               colNames = TRUE, rowNames = TRUE, sheetName = myContr[i],
               keepNA = TRUE, firstRow = TRUE) # Freezes the first row!
  }
}

summary.RP = rbind(colSums(results.RP == -1), # Cast to matrix
                   colSums(results.RP == 0),
                   colSums(results.RP == 1))
rownames(summary.RP) = c("Down", "NotSig", "Up")
summary.RP





# DEG Comparison -------------------------------------------------------------------------------------------------------
# Venn diagrams of DEGs from Limma and RP methods

# To suppress 'venn.diagram()' logging
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagrams
for (i in 1:length(myContr)) {
  for (j in c(1, -1)) {
    
    DEG.id.limma = rownames(results.limma)[which(results.limma[,i] == j)]
    DEG.id.RP = rownames(results.RP)[which(results.RP[,i] == j)]
    
    if (j == 1) {
      venn.sub = "UP-regulated DEGs"
    } else {
      venn.sub = "DOWN-regulated DEGs"
    }
    
    if (length(DEG.id.limma) == 0 & length(DEG.id.RP) == 0) {
      cat("\nBoth the sets are empty in the contrast: ", myContr[i], " (", venn.sub, ")\n\n", sep = "")
      next # Skip the current iteration of the for loop without terminating it
    }
    
    venn.plot = venn.diagram(x = list(DEG.id.limma, DEG.id.RP),
                             filename = NULL, # to print just on screen
                             force.unique = TRUE,
                             main = myContr[i], main.cex = 2, main.fontface = "bold", main.fontfamily = "sans", # Title
                             sub = venn.sub, sub.fontfamily = "sans", # Subtitle
                             lwd = 2, lty = "blank", fill = myColors[1:2], # circles
                             cex = 2, fontface = "bold", fontfamily = "sans", # numbers
                             category.names = c("Limma", "Rank Product"), # names
                             cat.cex = 2,
                             cat.fontface = "bold",
                             cat.default.pos = "outer",
                             cat.pos = c(-150, 150),
                             cat.dist = c(0.055, 0.055),
                             cat.fontfamily = "sans")
    
    # Create a new canvas and draw the Venn
    grid.newpage()
    grid.draw(venn.plot)
    printPlots(paste("11 - Comparison Venn ", myContr[i], "_", strsplit(venn.sub, split = "-")[[1]][1], sep = ""))
  }
}




