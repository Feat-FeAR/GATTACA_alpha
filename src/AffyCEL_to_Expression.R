# Header Info ------------------------------------------------------------------
#
# Affymetrix CEL files --to--> expression matrix
# For Affymetrix 3'IVT and Gene/Exon ST Arrays
#
# a FeAR R-script - 28-Dec-2021
# 
# Based on:
# [1] "Homer MiniTutorial - Processing Affymetrix Gene Expression Arrays"
#      by Chris Benner et al.
#      http://homer.ucsd.edu/homer/
#      http://homer.ucsd.edu/homer/basicTutorial/affymetrix.html
#
# [2] "An end to end workflow for differential gene expression using Affymetrix
#      microarrays"
#      by Bernd Klaus and Stefanie Reisenauer
#      https://bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html
#
# NOTE
# On both Gene 1.0 ST and Exon 1.0 ST Arrays a probe set is more or less an exon.
# Doing the analysis on probe set level means to analyze signals for each exon.
# On the contrary, a 'transcript cluster' contains all the probe sets of a gene,
# and therefore can be used to measure gene expressions. Accordingly, to perform
# analyses on gene level, always use target="core" as summarization target in
# the following functions from oligo package and "Transcript Cluster Annotations"
# packages that contain all annotations of the gene each probe set belongs to.
#
# General Script for CEL-file normalization through RMA algorithm
#
#   - CEL file loading
#   - RMA normalization
#       1 - Background correcting (normexp)
#       2 - Normalizing (Quantile-Quantile)
#       3 - Calculating expression (ProbeSet Summarization)
#   - Remove invalid probes
#   - Add annotation
#   - Save expression matrix
#
# Symbol Legend in Document Outline:
#   (NONE)  This section can be skipped (from a syntactical point of view)
#   *       This section is mandatory for subsequent blocks
#   **      This section is mandatory and requires a manual editing by the user
#   ~       This section is mandatory, but alternative to the previous one
#





# * Set Paths ------------------------------------------------------------------
# Set i/o paths

# Desktop local path
myDesktop = paste(Sys.getenv("USERPROFILE"), "Desktop", sep = .Platform$file.sep)

# Input and output folder
# NOTE: This way of selecting a file only works within RStudio !!
CEL.Dir = rstudioapi::selectDirectory(caption = "Select .CEL-containing folder",
                                      label = "Select",
                                      path = myDesktop)
out.Dir = rstudioapi::selectDirectory(caption = "Select output folder",
                                      label = "Select",
                                      path = paste0(CEL.Dir, "/.."))
setwd(out.Dir)





# ** Variable Definition -------------------------------------------------------

# Boolean flag to choose between remote and local database
use.remote = TRUE

platform = "hgu133a"
# Choose among the following:
# hgu133a                       Affymetrix Human Genome U133 A Set
# hgu133b                       Affymetrix Human Genome U133 B Set
# hgu133plus2                   Affymetrix Human Genome U133 Plus 2.0 Array
# hugene10stprobeset            Affymetrix GeneChip Human Gene 1.0 ST Array
#                                 ProbeSet Annotations -> Exon Level
# hugene10sttranscriptcluster   Affymetrix GeneChip Human Gene 1.0 ST Array
#                                 Transcript cluster Annotations-> Gene Level

# Check Platform and set the exon.probes flag
if (platform %in% c("hgu133a", "hgu133b", "hgu133plus2")) {
  exon.probes = FALSE # Affy 3'IVT arrays
} else if (platform %in% c("hugene10stprobeset", "hugene10sttranscriptcluster")) {
  exon.probes = TRUE  # Affy ST arrays
} else {
  cat("\n")
  stop("Invalid or unsupported platform\n\n")
}

# Global options suitable for STALKER_Functions
# Erase them by typing: options(scriptName = NULL, ...)
options(scriptName = "CEL_2_Exprs",
        save.PNG.plot = TRUE,
        save.PDF.plot = FALSE,
        append.annot = FALSE)





# * Load Packages --------------------------------------------------------------

# Load Annotation Package
if (getOption("append.annot") && use.remote) {
  annot.db = paste0(platform, ".db")
  # library() converts its argument into a string unless you specify
  # the option character.only = TRUE
  library(annot.db, character.only = TRUE)
  cat("Loaded annotation: ", annot.db, "\n\n", sep = "")
}

library(oligo)  # For GeneChip/Exon Affymetrix Arrays
                # AND most popular 3'IVT Affymetrix platforms
# As a rule for 3'IVT Affymetrix Arrays: try oligo package first, then if it
# doesn't work, decomment the following line and try affy package.
#library(affy)  # Most suitable for older 3'IVT (e.g. Mouse 430, Human U133, ...)

library(affycoretools)  # Utility functions for Affymetrix GeneChips

# Collection of custom functions
# NOTE: This way of sourcing only works within RStudio !!
GATTACA.dir = dirname(rstudioapi::getSourceEditorContext()$path) 
source(file.path(GATTACA.dir, "STALKER_Functions.R", fsep = .Platform$file.sep))





# * Load CEL files -------------------------------------------------------------

celFiles = list.celfiles(CEL.Dir, full.name = TRUE)
as.matrix(basename(celFiles))

# Parser to CEL files
# This should also install-and-load the platform design package (e.g. pd.hg.u133a)
affyRaw = read.celfiles(celFiles) # ExpressionFeatureSet data structure

# Before-normalization quality control Boxplots
set.seed(101) # For reproducibility, because oligo::boxplot methods use a sample
              # of the probesets to produce the plot
if (exon.probes) {
  # For Gene ST and Exon ST arrays, oligo::boxplot() needs the 'target' option,
  # when used with a FeatureSet-like object. 'target' argument can take one of
  # the following values: "probeset", "core" (Gene/Exon), "full" (Exon),
  # "extended" (Exon), describing the summarization target.
  oligo::boxplot(affyRaw, target = "core", names = c(1:length(celFiles)),
                 xlab = "Sample Index", main = "Unnormalized Expression",
                 ylab = "log2(intensity)", transfo = log2) # Raw data to be logged
} else {
  oligo::boxplot(affyRaw, names = c(1:length(celFiles)), xlab = "Sample Index",
                 main = "Unnormalized Expression", ylab = "log2(intensity)",
                 transfo = log2) # Raw data to be logged
}
printPlots("Boxplots - Pre-Norm")

# Before-normalization quality control MA-plots
n = 4                 # Number of samples to MA-plot
m = length(celFiles)  # Total number of samples (i.e. expression matrix columns)
# Sample subset by sampling n (out of m) chips
if (n >= m) {
  sample.index = round(seq(from = 1, to = m, length = m))
} else {
  sample.index = round(seq(from = 1, to = m, length = n))
}
# MAplot function returns 1 chip vs the pseudo-median reference chip
for (i in sample.index) {
  oligo::MAplot(affyRaw, which = i, pairs = FALSE)
  printPlots(paste0("MA-Plot_", i, " - Pre-Norm"))
}





# * RMA Normalization ----------------------------------------------------------

if (exon.probes) {
  # The parameter 'target' (only for Gene ST and Exon ST) defines the degree of
  # summarization: "core" is the default option which use transcript clusters
  # containing "safely" annotated genes. For summarization on the exon level
  # (not recommended for Gene arrays), use "probeset".
  eset = rma(affyRaw, target = "core") # ExpressionSet data structure
} else {
  eset = rma(affyRaw) # ExpressionSet data structure
}

expressionMatrix = exprs(eset)
d = show.data(expressionMatrix)

# After-normalization quality control Boxplots
oligo::boxplot(eset, names = c(1:length(celFiles)), xlab = "Sample Index",
               main = "RMA-Normalized Expression", ylab = "log2(intensity)",
               transfo = identity) # Already logged (by RMA)
printPlots("Boxplots - Post-Norm")

# After-normalization quality control MA-plots
# Single samples vs pseudo-median reference chip
for (i in sample.index) {
  oligo::MAplot(eset, which = i, pairs = FALSE)
  printPlots(paste0("MA-Plot_", i, " - Post-Norm"))
}





# Remove Invalid Probes --------------------------------------------------------

# Remove Affymetrix control probes
probes.before = dim(eset)[1]
if (exon.probes) {
  eset = getMainProbes(eset, level = "core")
} else {
  # ^ anchor to match the start of string (see regular expressions)
  ctrl.index = grep("^AFFX", rownames(eset))
  # To prevent an integer(0) from deleting all probes
  if (length(ctrl.index) > 0) {
    eset = eset[-ctrl.index,]
  }
}
probes.after = dim(eset)[1]
discarded = probes.before - probes.after
cat("\n", discarded, " Affymetrix control probes (",
    round(1e2*discarded/probes.before, digits = 2),
    "%) have been detected and discarded\n\n", sep = "")

# Check for missing values (NA) and NaN entries
expressionMatrix = exprs(eset)
any(is.na(expressionMatrix))
any(is.nan(expressionMatrix))





# Add Annotations --------------------------------------------------------------

if (getOption("append.annot")) {
  annot = create.annot(platform, remote = use.remote)
} else {
  annot = NULL
  cat("\nNo annotation loaded\n\n", sep = "")
}

# Count missing annotations in the complete db
if (getOption("append.annot")) {
  
  if (use.remote) {naSymb = "NA"} else {naSymb = "---"}
  
  nc = dim(annot)[2]
  notMap = matrix(0, nrow = 2, ncol = nc,
                  dimnames = list(c("Not Mapped","%"), colnames(annot)))
  for (i in 1:nc) {
    notMap[1,i] = length(which(annot[,i] == naSymb))
    notMap[2,i] = round(notMap[1,i]/dim(annot)[1]*1e2, digits = 2)
  }
  notMap
}

myFrame = data.frame(expressionMatrix) # Matrix to data frame
myFrame = appendAnnotation(myFrame, annot)
show.data(myFrame)

# Count missing annotations in the actual expression matrix
if (getOption("append.annot")) {
  for (i in 1:nc) {
    notMap[1,i] = length(which(myFrame[,i] == naSymb))
    notMap[2,i] = round(notMap[1,i]/dim(myFrame)[1]*1e2, digits = 2)
  }
  notMap
}





# Save Expression Matrix -------------------------------------------------------

# Without annotations
write.csv(expressionMatrix, row.names = TRUE,
          file = paste0("RMA-normalized_logExpression_", platform, ".csv"))
cat("\nExpression matrix--without annotations--has been saved as CSV file in\n",
    out.Dir, "\n\n", sep = "")

# With annotations
if (getOption("append.annot")) {
  write.csv(myFrame, row.names = TRUE,
            file = paste0("RMA-normalized_logExpression_", platform, "_annotated.csv"))
  cat("\nExpression matrix--with annotations--has been saved as CSV file in\n",
      out.Dir, "\n\n", sep = "")
} else {
  cat("\nNo annotation loaded",
      "\nExpression matrix with annotations cannot be saved\n\n", sep = "")
}



