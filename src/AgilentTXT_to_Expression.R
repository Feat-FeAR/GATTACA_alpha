# Header Info ------------------------------------------------------------------
#
# Agilent TXT raw data --to--> expression matrix
#
# a FeAR R-script - 27-Jan-2023
#
# NOTE
# read.maimages() function from limma requires as its first argument a data frame
# containing (at least) a column called 'FileName' with the names of the raw-data
# files to be used for the analysis. For Agilent arrays such a list is usually
# given in the form of a tab-delimited file named "Targets.txt", possibly
# containing other kind of information about the experimental design. However,
# this script ignores any Targets.txt file found within the input directory and
# builds its own target list run-time.
#
# Script outline:
#   - Raw data loading
#   - 'normexp' background correction
#   - Quantile-Quantile interarray normalization
#   - Negative control probe evaluation
#   - Invalid probes removal
#   - Add annotation
#   - Expression matrix saving
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
myDesktop <- paste(Sys.getenv("USERPROFILE"), "Desktop",
                   sep = .Platform$file.sep)

# Input and output folder
# NOTE: This way of selecting a file only works within RStudio !!
TXT.Dir <- rstudioapi::selectDirectory(caption = "Select TXT-containing folder",
                                      label = "Select",
                                      path = myDesktop)
out.Dir <- rstudioapi::selectDirectory(caption = "Select output folder",
                                      label = "Select",
                                      path = paste0(TXT.Dir, "/.."))
setwd(out.Dir)





# ** Variable Definition -------------------------------------------------------

# Boolean flag to choose between remote and local database
use.remote <- TRUE

# Number of MA-plots to be drawn for each Agilent.QC.Plots
maplot.num <- 8

# Global options suitable for STALKER_Functions
# Erase them by typing: options(scriptName = NULL, ...)
options(scriptName = "TXT_2_Exprs",
        save.PNG.plot = TRUE,
        save.PDF.plot = FALSE,
        append.annot = FALSE)





# * Load Packages --------------------------------------------------------------

# Load Annotation Package
if (getOption("append.annot") && use.remote) {
  annot.db <- paste0(platform, ".db")
  # library() converts its argument into a string unless you specify
  # the option character.only = TRUE
  library(annot.db, character.only = TRUE)
  cat("Loaded annotation: ", annot.db, "\n\n", sep = "")
}

library(limma)      # Empirical Bayes Method for Differential Expression
                    # here used for low-level analysis only

# Collection of custom functions
# NOTE: This way of sourcing only works within RStudio !!
GATTACA.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(file.path(GATTACA.dir, "STALKER_Functions.R", fsep = .Platform$file.sep))

# Array platform
platform <- array_platform_selector("Agilent")





# * Load TXT raw data ----------------------------------------------------------

# List all .TXT files in TXT.Dir... except possible "Targets.txt" file
txtFiles <- list.files(path = TXT.Dir, full.names = TRUE,
                      pattern = "*.(txt|TXT)")
target.index <- grep("targets.txt", txtFiles, ignore.case = TRUE)
if (length(target.index)>0) {txtFiles <- txtFiles[-target.index]}
as.matrix(basename(txtFiles))

# Build a Target File run-time
targets <- data.frame(SampleNumber = c(1:length(txtFiles)), FileName = txtFiles)
targets

# Convert the data to an EListRaw object (data object for single-channel arrays)
# Specify green.only=TRUE for one-channel data and retain information about
# background via gIsWellAboveBG
raw <- read.maimages(targets, source = 'agilent.median', green.only = TRUE,
                     other.columns = 'gIsWellAboveBG')

# Before-normalization quality control plots (raw data to be log2-transformed)
Agilent.QC.Plots(log2(as.matrix(raw)), targets, "Pre-Norm",
                 maplots = maplot.num)





# Sample selection -------------------------------------------------------------

# Exclude corrupted, invalid, or blank arrays
targets <- targets[-4,]
raw <- read.maimages(targets, source = 'agilent.median', green.only = TRUE,
                    other.columns = 'gIsWellAboveBG')





# * Normalization --------------------------------------------------------------
# Background subtraction and inter-array normalization

# Subtract the background signal
# Offset to shrunk towards zero log-ratios at the lower intensities
# see 'anti-fanning in MA-Plot' for more info on this.
ofst <- 50
raw_BGcorrected <- backgroundCorrect(raw, method = "normexp", offset = ofst)

# BG-subtracted quality control plots (Data still to be log2-transformed)
Agilent.QC.Plots(log2(as.matrix(raw_BGcorrected)), targets, "BG-subtracted",
                 maplots = maplot.num)

# Quantile-normalize and log-transform the data
raw_BGandNormalized <- normalizeBetweenArrays(raw_BGcorrected,
                                              method = "quantile")

# After-normalization quality control plots (Already log2-transformed data)
Agilent.QC.Plots(as.matrix(raw_BGandNormalized), targets, "Post-Norm",
                 maplots = maplot.num)

# Extract the expression matrix from EListRaw object and rename columns
expressionMatrix <- raw_BGandNormalized$E
colnames(expressionMatrix) <- paste0("Sample_", targets$SampleNumber)

# Show preview
d <- show_data(expressionMatrix)





# Negative Control -------------------------------------------------------------
# Check the negative control probe (-)3xSLv1 to estimate the log2_intensity of
# unhybridized spots, to be used as the threshold value for filtering in GATTACA

neg.ctrl <- raw_BGandNormalized$genes$ControlType == -1
neg.id <- unique(raw_BGandNormalized$genes$ProbeName[neg.ctrl])
#expressionMatrix[neg.ctrl,]

cat("\n", sum(neg.ctrl), "Negative-Control probes have been found,\n",
    "corresponding to", length(neg.id), "unique probe(s):", neg.id, "\n",
    "Mean value from", length(expressionMatrix[neg.ctrl,]),
    "unhybridized spots:", mean(expressionMatrix[neg.ctrl,]), "\n\n", sep = " ")





# Remove Invalid Probes --------------------------------------------------------

# Filter out control probes
# NOTE: abs() is used here to include also the NegativeControl probe (-)3xSLv1
ctrlProbes <- abs(raw_BGandNormalized$genes$ControlType) == 1
cat("\n", sum(ctrlProbes), " control probes detected\n\n", sep = "")
# Get names of Control Probes detected
unique(raw_BGandNormalized$genes$ProbeName[ctrlProbes])
#unique(as.data.frame(raw_BGandNormalized)[ctrlProbes, 7]) # Alternatively

# Filter out those probes that fail the 'gIsWellAboveBG' check
isExpr <- rep(TRUE, dim(raw_BGandNormalized)[1])
if (FALSE) { # No second filter - Default
  # Well above BG in at least 50% of the samples
  isExpr <- rowSums(raw_BGandNormalized$other$gIsWellAboveBG > 0) >= d[2]/2
  sum(isExpr)
}

raw_BGandNormalized_filt <- raw_BGandNormalized[!ctrlProbes & isExpr, ]

# Probe Summarization: replace values of replicate probes with their mean
# Probe_IDs are used to identify the replicates and then assigned to row names
raw_BGandNormalized_filt_mean <- avereps(raw_BGandNormalized_filt,
                                        ID = raw_BGandNormalized_filt$genes$ProbeName)
pre <- dim(raw_BGandNormalized_filt)
post <- dim(raw_BGandNormalized_filt_mean)
cat("\n", (pre - post)[1], " replicate probes have been averaged\n\n", sep = "")
expressionMatrix <- raw_BGandNormalized_filt_mean$E
colnames(expressionMatrix) <- paste0("Sample_", targets$SampleNumber)
d <- show_data(expressionMatrix)





# Add Annotations --------------------------------------------------------------

if (getOption("append.annot")) {
  annot <- create.annot(platform, remote = use.remote)
} else {
  annot <- NULL
  cat("\nNo annotation loaded\n\n", sep = "")
}

# Count missing annotations in the complete db
if (getOption("append.annot")) {
  
  if (use.remote) {naSymb = "NA"} else {naSymb = "---"}
  
  nc <- dim(annot)[2]
  notMap <- matrix(0, nrow = 2, ncol = nc,
                  dimnames = list(c("Not Mapped","%"), colnames(annot)))
  for (i in 1:nc) {
    notMap[1,i] <- length(which(annot[,i] == naSymb))
    notMap[2,i] <- round(notMap[1,i]/dim(annot)[1]*1e2, digits = 2)
  }
  notMap
}

myFrame <- data.frame(expressionMatrix) # Matrix to data frame
myFrame <- appendAnnotation(myFrame, annot)
# Show preview
d <- show_data(myFrame)

# Count missing annotations in the actual expression matrix
if (getOption("append.annot")) {
  for (i in 1:nc) {
    notMap[1,i] <- length(which(myFrame[,i] == naSymb))
    notMap[2,i] <- round(notMap[1,i]/dim(myFrame)[1]*1e2, digits = 2)
  }
  notMap
}





# Save Expression Matrix -------------------------------------------------------

# Without annotations
write.csv(expressionMatrix, row.names = TRUE,
          file = paste0("Normalized_logExpression_", platform, ".csv"))
cat("\nExpression matrix--without annotations--has been saved as CSV file in\n",
    out.Dir, "\n\n", sep = "")

# With annotations
if (getOption("append.annot")) {
  write.csv(myFrame, row.names = TRUE,
            file = paste0("Normalized_logExpression_", platform, "_annotated.csv"))
  cat("\nExpression matrix--with annotations--has been saved as CSV file in\n",
      out.Dir, "\n\n", sep = "")
} else {
  cat("\nNo annotation loaded",
      "\nExpression matrix with annotations cannot be saved\n\n", sep = "")
}



