# Header Info ------------------------------------------------------------------
# STALKER Functions
#
# A Collection of R functions to be used with:
#     GATTACA_alpha
#
# a //FeAR// R-script - 25-Jan-2023
#

#'------------------------------------------------------------------------------
#' @description A custom version of the classical `head()` that prints the upper
#'              leftmost corner of a data set, also showing row names and
#'              controlling for possible out-of-bounds exceptions. Compared to
#'              `head()`, `show_data()` *always displays data by columns*, even
#'              in the case of vectors (i.e., one-dimensional arrays). Also,
#'              `show_data` prints and returns the dimensions of the data set,
#'              together with a custom heading label.
#' 
#' @param dataset Data frame, matrix, or vector to print.
#' @param name Explanatory name to print in the heading (useful when logging).
#' @param rows Maximum number of rows to display.
#' @param cols Maximum number of columns to display.
#' 
#' @returns A vector containing `dataset` dimensions.
#' 
#' @author //FeAR//
#'------------------------------------------------------------------------------
show_data <- function(dataset, name = NULL, rows = 10, cols = 5)
{
  # Whenever possible, 'dataset' is converted to a data frame to be passed to
  # dim() function and displayed column-wise preserving the type specificity of
  # each column. Notably, duplicated names are allowed in both named vectors and
  # matrices, but not in data frames. So, when converting to data frame, names
  # of named vectors are discarded (or better "reassigned to integers) if and
  # only if duplicates are present, while possible duplicated row names from
  # matrices are disambiguated by progressively appending the suffixes '.1',
  # '.2', '.3'... to the original row names!! For all these reasons, matrices
  # and named vectors with duplicated names will be converted to matrices
  # instead of data frames before being printed, in order to always preserve
  # original row names. Note that all these conversions only concern the printed
  # copy of 'dataset', while the type of the original data is unaffected.
  if ((is.vector(dataset) & sum(duplicated(names(dataset))) > 0) |
      (is.matrix(dataset) & sum(duplicated(row.names(dataset))) > 0)) { 
    # If 'dataset' is a named vector with duplicated names or a matrix with
    # duplicated row names, convert it to a matrix...
    dataset <- as.matrix(dataset)
  } else {
    # ...otherwise convert it to a data frame.
    dataset <- as.data.frame(dataset)
  }
  d <- dim(dataset)
  cat("\nDataset", name, "dimensions:", d[1], "x", d[2], "\n\n")
  
  rows <- min(d[1], rows)
  cols <- min(d[2], cols)
  
  # 'print' because automatic printing is turned off in loops (and functions)
  # NOTE: if you only return a data frame subset of one column, R will drop the
  # names by default. To avoid this use the drop=FALSE option.
  print(dataset[1:rows, 1:cols, drop = FALSE])
  
  return(d)
}

#'------------------------------------------------------------------------------
#' @description Hypergeometric test function.
#' 
#' @param k Number of hits in the experimental set.  --OR--> Intersection
#' @param n Size of the experimental set.            --OR--> Set A
#' @param K Number of possible hits in the universe. --OR--> Set B
#' @param N Size of the universe.                    --OR--> Background
#'
#' @returns A 1-by-5 data frame containing enrichment statistics. Many
#'          single-row data frames can be easily stacked using `rbind()`.
#'          
#' @author //FeAR//
#'------------------------------------------------------------------------------
hgt <- function(k, n, K, N = 1e4)
{
  pval <- phyper(k-1, K, N-K, n, lower.tail = FALSE) # p-value
  expect <- n*(K/N) # Expected value
  stdev <- sqrt(n*(K/N)*((N-K)/N)*((N-n)/(N-1))) # Standard Deviation
  FE <- k/(n*(K/N)) # Fold Enrichment
  
  stats <- data.frame(Hits = k,
                      Expected_Value = expect,
                      SD = stdev,
                      Fold_Enrichment = FE,
                      p.value = pval)
  return(stats)
}

#'-NEW--------------------------------------------------------------------------
#' @description A minimal graphical interface to retrieve the name of a suitable
#'              annotation database (a `.db` R-package or a `GPL` GEO platform
#'              record) starting from the selection of the name of a microarray
#'              chip model.
#'
#' @param filt A string used to subset the list of platform annotations among
#'             which to choose (i.e., a `grep` filter working on array names).
#'             It can target array manufacturer ("Affymetrix", "Agilent"),
#'             or the source of annotation ("Biocondutor", "GPL" for GEO).
#'
#' @returns The name of the database corresponding to the platform chosen by the
#'          user (to be used with `create.annot()` function).
#'
#' @example array_platform_selector() |> array_create_annot(platform) |> head()
#'           
#' @author //FeAR//
#'------------------------------------------------------------------------------
array_platform_selector <- function(filt = "All")
{
  # Name of the Bioconductor db package or GPL GEO platform record
  db_BCorGPL <- c("hgu133a",
                  "hgu133b",
                  "hgu133plus2",
                  "hugene10stprobeset",
                  "hugene10sttranscriptcluster",
                  "",
                  "hgug4112a",
                  "GPL6480",
                  "HsAgilentDesign026652",
                  "GPL22763",
                  "GPL19072")
  
  # Array full-length name or `title` value from @header slot of a GPL object
  # NOTE_1: GEOquery::getGEO() function needs FTP
  # NOTE_2: Only Agilent 'Probe Name Versions' (NO 'Gene Symbol Versions') GPLs
  #         have been included in this list. 
  long_names <- c("Affymetrix Human Genome U133 A Set - from Bioconductor",
                  "Affymetrix Human Genome U133 B Set - from Bioconductor",
                  "Affymetrix Human Genome U133 Plus 2.0 Array - from Bioconductor",
                  "Affymetrix GeneChip Human Gene 1.0 ST Array - Exon Level - from Bioconductor",
                  "Affymetrix GeneChip Human Gene 1.0 ST Array - Gene Level - from Bioconductor",
                  "---",
                  "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F - from Bioconductor",
                  "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F - from GPL6480",
                  "Agilent-026652 Whole Human Genome Microarray 4x44K v2 G4845A - from Bioconductor",
                  "Agilent-039714 LincRNA SurePrint G3 Human GE 8x60K Microarray PVD 028004 - from GPL22763",
                  "Agilent-052909 CBC_lncRNAmRNA_V3 - from GPL19072")
  
  if (filt != "All") {
    filt_index <- grep(filt, long_names, ignore.case = TRUE)
    long_names_sub <- long_names[filt_index]
    db_BCorGPL_sub <- db_BCorGPL[filt_index]
  } else {
    long_names_sub <- long_names
    db_BCorGPL_sub <- db_BCorGPL
  }
  
  platform_index <- menu(long_names_sub, title = "Choose platform annotation",
                         graphics = TRUE)
  
  if (platform_index == 0) {
    cat("\nNo platform selected!\n\n")
  } else {
    cat(paste("\nSelected platform:\n", long_names_sub[platform_index], "\n\n"))
    
    return(db_BCorGPL_sub[platform_index])
  }
}

#'-NEW--------------------------------------------------------------------------
#' @description Create the annotation data frame starting from a database name
#'              as returned by the `array_platform_selector()`. This function
#'              implements a minimal graphical interface allowing the user to
#'              select the number and the type of features to be used as
#'              annotation.
#' 
#' @param platform Affymetrix/Agilent platform annotation database.
#' @param collapsing Boolean flag to choose whether to collapse by unique Probe
#'                   ID in the case of annotation packages from Bioconductor
#'                   (GPL records from GEO are already collapsed).
#' 
#' @returns A data frame containing for each probe of the platform an number of
#'          features selected by the user.
#'
#' @example annot <- array_create_annot("hgu133plus2", collapsing = TRUE)
#'          missing_report(annot)
#'          d <- show_data(annot)
#'          
#' @author //FeAR//
#'------------------------------------------------------------------------------
array_create_annot <- function(platform, collapsing = FALSE)
{
  # Download the annotation matrix from GEO and subset by columns of interest  
  if (grepl("GPL", platform)) {
    
    # Download annotation, extract the `dataTable`, print some information
    GEO_GPL <- GEOquery::getGEO(platform)
    annot_full <- GEOquery::Table(GEO_GPL)
    cat("\nLoaded annotation: ", platform,
        " (Last Update: ", GEO_GPL@header$last_update_date, ")", sep = "")
    
    # Retrieve all the columns of the database for subsequent feature selection
    cols <- colnames(annot_full)
    feats <- svDialogs::dlg_list(choice = cols,
                                 preselect = c("ENSEMBL_ID",  # Ensembl
                                               "GENE",        # NCBI Entrez Gene
                                               "GENE_SYMBOL", # HGNC Gene Symbol
                                               "GENE_NAME"),  # Gene Name
                                 multiple = TRUE,
                                 title = "Select multiple features")$res
    
    # Always use Probe IDs as keys (GPL matrices are already collapsed by unique
    # Probe IDs)
    ids <- annot_full[,1] # "Probe ID" feature should always be the first column
    cat("\n", length(ids), " unique Probe_IDs retrieved from ",
        platform," GEO platform record\n\n", sep = "")
    
    # Subset (including Probe IDs) and return
    annot <- annot_full[,c(cols[1],feats)]
    return(annot)
    
  } else {
    
    # Load Annotation Database from Bioconductor and retrieve columns of interest
    annot_db <- paste0(platform, ".db")
    # library() converts its argument into a string unless you specify the
    # option character.only = TRUE
    library(annot_db, character.only = TRUE)
    library(AnnotationDbi)
    # Print some information
    cat("\nLoaded annotation: ", annot_db,
        " (ver.: ", toString(packageVersion(annot_db)), ")",
        " [date: ", toString(packageDate(annot_db)), "]", sep = "")
    
    # The evaluated expression of the (unquoted) annot_db (e.g., hgu133a.db)
    evaluating_db <- eval(parse(text = annot_db))
    
    # Retrieve all the columns of the data base, except "PROBEID" that will be
    # used as key and so it must be unique
    cols <- columns(evaluating_db)
    cols <- cols[! cols %in% "PROBEID"]
    
    # Feature selection
    feats <- svDialogs::dlg_list(choice = cols,
                                 preselect = c("ENSEMBL",
                                               "ENTREZID",
                                               "SYMBOL",
                                               "GENENAME"),
                                 multiple = TRUE,
                                 title = "Select multiple features")$res
    
    # Always use Probe IDs as keys (although they cannot be row names since, in
    # general, they will not be unique after feature retrieval)
    ids <- keys(evaluating_db, keytype = "PROBEID")
    cat("\n", length(ids), " unique Probe_IDs retrieved from ",
        annot_db,"\n", sep = "")
    annot_long <- select(evaluating_db,
                         keys = ids,
                         columns = feats,
                         keytype = "PROBEID")
    cat("1:many mapping resulted in a ",
        dim(annot_long)[1], " x ", dim(annot_long)[2],
        " annotation data frame", sep = "")
    
    if (collapsing) {
      del = " /// " # Affymetrix-style delimiter
      # Always use Probe IDs as keys
      annot_collapse <- data.frame(PROBEID = keys(evaluating_db,
                                                  keytype = "PROBEID"))
      for (feat in feats) {
        annot_collapse <- merge(annot_collapse,
                                aggregate(annot_long[feat],
                                          by = annot_long["PROBEID"],
                                          function(...){paste(unique(...),
                                                              collapse = del)}),
                                by.x = "PROBEID", by.y = "PROBEID", all.x = TRUE)
      }
      cat("\n...now reduced to ", dim(annot_collapse)[1], " x ",
          dim(annot_collapse)[2], " after collapsing on Probe_IDs\n\n", sep = "")
      return(annot_collapse)
      
    } else {
      cat("\n\n")
      return(annot_long)
    }
  }
}

#'-NEW--------------------------------------------------------------------------
#' @description This function takes a data frame, searches its columns for many
#'              common missing-value placeholders, and finally prints a report
#'              report of the number and the relative amount of missing values
#'              detected column-wise. In particular, it searches for `NA`,
#'              `"NA"`, and patterns of whitespaces (`\s`), hyphens (`-`),
#'              slashes (`/`), including empty fields and a possible
#'              user-defined sequence.
#' 
#' @param dataFrame Data frame or matrix to be scanned for NAs.
#' @param naSymb A string containing the user-defined sequence for NAs.
#'
#' @example array_platform_selector() |> array_create_annot() |> missing_report()
#' 
#' @author //FeAR//
#'------------------------------------------------------------------------------
missing_report <- function(dataFrame, naSymb = "")
{
  cols <- colnames(dataFrame)
  missing_data <- matrix(0, nrow = 2, ncol = length(cols),
                         dimnames = list(c("Not Mapped","%"), cols))
  
  # Search for NAs, "NA"s and patterns of whitespaces (\s), hyphens (-),
  # slashes (/), including empty fields (*) and a user-defined sequence
  for (feat in cols) {
    missing_data[1,feat] <- sum(
      is.na(dataFrame[,feat]) |
        dataFrame[,feat] == "NA" |
        grepl("^(\\s|-|/)*$", dataFrame[,feat]) |
        dataFrame[,feat] == naSymb)
    missing_data[2,feat] <- round((missing_data[1,feat]/dim(dataFrame)[1])*1e2,
                                  digits = 2)
  }
  cat("\n")
  print(missing_data)
  cat("\n")
}

#'------------------------------------------------------------------------------
#' @description Save a graphical output to '<folderPrefix> Figures'
#'              sub-directory.
#'
#' @param figureName Name of the output file (without extension).
#' @param folderPrefix Prefix for naming the saving subfolder (default to the
#'                     name of the parent script sourcing this function).
#' @param PNG Boolean: T to print the currently displayed figure in PNG format.
#' @param PDF Boolean: T to print the currently displayed figure in PDF format.
#'
#' @author //FeAR//
#'------------------------------------------------------------------------------
printPlots = function(figureName, folderPrefix = getOption("scriptName"),
                      PNG = getOption("save.PNG.plot"),
                      PDF = getOption("save.PDF.plot"))
{
  flag = FALSE # A dummy flag to insert a couple of 'new lines' in case of WARNINGs
  
  # Check argument values
  # NOTE: considering that getOption("...") returns NULL for undefined arguments,
  #       IFs are evaluated only when:
  #       the corresponding global option is not defined
  #         AND
  #       no argument is passed runtime
  if (is.null(PNG)) { 
    PNG = TRUE
    flag = TRUE
    cat("\nWARNING: \'save.PNG.plot\' option defaulted to TRUE")
  }
  if (is.null(PDF)) {
    PDF = TRUE
    flag = TRUE
    cat("\nWARNING: \'save.PDF.plot\' option defaulted to TRUE")
  }
  if (is.null(folderPrefix)) {
    figSubFolder = "Figures"
  } else {
    figSubFolder = paste(folderPrefix, "Figures", sep = " ")
  }
  
  fullName = file.path(figSubFolder, figureName, fsep = .Platform$file.sep)
  
  if (!file.exists(figSubFolder) && (PNG || PDF)) {
    dir.create(figSubFolder)
    flag = TRUE
    cat("\nNew folder '", figSubFolder, "' has been created in the current WD",
        sep = "")
  }
  if (PNG) { # invisible(capture.output()) to suppress automatic output to console
    invisible(capture.output(
      dev.print(device = png, filename = paste0(fullName, ".png"),
                width = 820, height = 600)))
  }
  if (PDF) {
    invisible(capture.output(
      dev.print(device = pdf, paste0(fullName, ".pdf"))))
  }
  if (flag) {
    cat("\n\n")
  }
}

#'------------------------------------------------------------------------------
#' @description Create the annotation data frame.
#' 
#' @param platform Affymetrix/Agilent platform name (only used for remote db).
#' @param remote Boolean flag to choose between remote and local database.
#' 
#' @returns A data frame containing Accession Number, GeneSymbol, and
#'          Description for each probe of the platform.
#'
#' @author //FeAR//
#'------------------------------------------------------------------------------
create.annot <- function(platform, remote = TRUE)
{
  if (remote) { # Use remote annotation
    
    # Load Annotation Package
    annot_db <- paste0(platform, ".db")
    # library() converts its argument into a string unless you specify the
    # option character.only = TRUE
    library(annot_db, character.only = TRUE)
    
    # Decomment the line below for a list of the available remote annotations
    # and additional information about version and time stamps
    #eval(parse(text = paste0(platform, "()")))
    
    acc <- eval(parse(text = paste0(platform, "ACCNUM")))
    symb <- eval(parse(text = paste0(platform, "SYMBOL")))
    gname <- eval(parse(text = paste0(platform, "GENENAME")))
    # NOTE: In general, probe sets map to other fields as one-to-many. In order
    #       to keep Probe IDs unique when creating the data frame, we need first
    #       to make each element of the lists `contents(...)` a single string
    #       (applying 'paste' together with 'collapse' option).
    del <- " /// " # Affy-style delimiter
    annot <- data.frame(Accession = sapply(contents(acc), paste, collapse = del),
                        GeneSymbol = sapply(contents(symb), paste, collapse = del),
                        Description = sapply(contents(gname), paste, collapse = del))
    cat("\nLoaded annotation: ", annot_db, sep = "")
    
  } else { # Use local annotation
    
    library(openxlsx) # Read, Write, and Edit .xlsx (Excel) Files
    myDesktop <- paste(Sys.getenv("USERPROFILE"), "Desktop",
                      sep = .Platform$file.sep)
    
    # NOTE: This way of selecting a file path only works within RStudio !!
    db.file <- rstudioapi::selectFile(caption = "Select Annotation File",
                                     label = "Select",
                                     path = myDesktop,
                                     filter = "All Files (*)",
                                     existing = TRUE)
    # Read as Data Frame
    annot <- read.xlsx(db.file, colNames = TRUE, rowNames = TRUE, sep.names = "_")
    # Decomment below for a list of the available local annotations
    #as.matrix(colnames(annot))
    annot <- annot[,c("Representative_Public_ID", "Gene_Symbol", "Gene_Title")]
  }
  
  cat("\nA ", dim(annot)[1], " x ", dim(annot)[2],
      " annotation dataframe has been loaded\n\n", sep = "")
  
  return(annot)
}

#'------------------------------------------------------------------------------
#' @description ppend annotation to genes and sort (do nothing if
#'              do.the.job == FALSE).
#' 
#' @param gene.stat The table of genes, usually a DEG summary-statistic
#'                  top-table (or an expression matrix).
#' @param ann The matrix containing the annotation data.
#' @param do.the.job FALSE to skip the appending task by global settings,
#'                   without the need for an external IF.
#' @param sort.by The name or index of the column used to sort the final dataset
#'
#' @returns The annotated and sorted data frame passed as input.
#'
#' @author //FeAR//
#'------------------------------------------------------------------------------
appendAnnotation = function(gene.stat, ann,
                            do.the.job = getOption("append.annot"),
                            sort.by = 1)
{
  # Check argument values
  if (is.null(do.the.job)) {
    do.the.job = TRUE
    cat("\nWARNING: \'append.annot\' option defaulted to TRUE\n\n")
  }
  
  if (do.the.job) {
    
    # 'merge' function to merge two matrix-like objects horizontally and cast to
    # data frame (right outer join).
    # NOTE: both gene.stat and ann are supposed to have the Probe_IDs as rownames
    
    # To merge two data frames horizontally by one or more common key variables:
    #  - inner join (default): Return only the rows that have matching keys in both
    #     the tables (~ intersection)
    #  - outer join (all = T): Return all rows from both the tables, joining the
    #     records that have matching (~ union)
    #  - left outer (all.x = T): Return all rows from the left table, and any rows
    #     with matching keys from the right table
    #  - right outer (all.y = T): Return all rows from the right table, and any rows
    #     with matching keys from the left table
    #  - cross join (by = NULL): Return the Cartesian product
    
    joined = merge(ann, gene.stat,
                   by.x = "row.names", by.y = "row.names", all.y = TRUE)
    rownames(joined) = joined[,1]
    gene.stat = joined[,-1]
    
    # Re-sort the data frame by the content of 'sort.by' column ('sort.by' can
    # be either a number or a column name).
    gene.stat = gene.stat[order(gene.stat[,sort.by]),]
  }
  return(gene.stat)
}

#'------------------------------------------------------------------------------
#' @description Return basics descriptive statistics of a single gene, by group
#'              label.
#' 
#' @param gene Numeric vector or single-row data frame from gene expression
#'             matrix.
#' @param gr Group names.
#' @param des Experimental design (full design mode vector).
#' @param prec Decimal precision.
#'
#' @returns A data frame containing the statistics of interest for each gene of
#'          `gene`.
#'
#' @author //FeAR//
#'------------------------------------------------------------------------------
descStat1G = function(gene, gr, des, prec = 4)
{
  # Define a new empty data frame
  stat.frame = data.frame(GROUP = character(),
                          n = integer(),
                          MEAN = double(),
                          VAR = double(),
                          SD = double(),
                          SEM = double(),
                          stringsAsFactors = FALSE)
  
  for (i in 1:length(gr)) {
    
    n.gene = as.numeric(gene[des == i]) # Downcast to numeric vector
    
    stat.frame[i,1] = gr[i]
    stat.frame[i,2] = sum(des == i)
    stat.frame[i,3] = round(mean(n.gene), digits = prec)
    stat.frame[i,4] = round(var(n.gene), digits = prec)
    stat.frame[i,5] = round(sd(n.gene), digits = prec)
    stat.frame[i,6] = round(sd(n.gene)/sqrt(sum(des == i)), digits = prec) # SEM
  }
  return(stat.frame)
}

#'------------------------------------------------------------------------------
#' @description Plot single gene comparison chart.
#' 
#' @param exp.mat Expression matrix (as data frame).
#' @param gr Group names.
#' @param des Experimental design (full design mode vector).
#' @param gois Genes of interest by probe (char vector).
#' @param chart.type "BP" (Box Plot), "BC" (Bar Chart), or "MS" (Mean & SEM).
#' @param ann Optional annotation data frame.
#'
#' @author //FeAR//
#'------------------------------------------------------------------------------
singleGeneView = function(exp.mat, gr, des, gois, chart.type = "BP", ann = NULL)
{
  geo = switch(chart.type,
               "BP" = "point",
               "BC" = "bar",
               "MS" = "crossbar")
  
  for (i in 1:length(gois)) {
    
    var.expr = as.numeric(exp.mat[gois[i],]) # Downcast to vector
    var.groups = gr[des]
    sgex = data.frame(var.expr, var.groups) # Single Gene Expression Data Frame
    sgs = descStat1G(exp.mat[gois[i],], gr, des, 6) # Single Gene Summary Data Frame
    
    if (is.null(ann)) {
      gene.symb = ""
    } else {
      gene.symb = paste(ann[gois[i], grepl("Symbol", colnames(ann))], " - ", sep = "")
    }
    
    if (chart.type == "BP") {
      
      print( # NOTICE: When in a for loop, you have to explicitly print your resulting ggplot object
        ggplot(data = sgex, aes(var.groups, var.expr)) +
          theme_bw(base_size = 15, base_rect_size = 1.5) +
          xlab("Group") + # In the following functions, when data=NULL (default), the data is inherited from ggplot()
          ylab("log2 Expression") +
          ggtitle(label = "Box Plot with Jitter", subtitle = paste(gene.symb, "Probe ID: ", gois[i], sep = "")) +
          geom_boxplot(width = 0.5, size = 0.5, notch = TRUE, outlier.shape = NA) +
          stat_summary(fun = "mean", geom = geo, color = "red3", size = 2) +
          geom_jitter(position = position_jitter(width = 0.1, height = 0, seed = 123), size = 1.5))
      
    } else if (chart.type == "BC" | chart.type == "MS") {
      
      print(
        ggplot(data = sgex, aes(var.groups, var.expr)) +
          theme_bw(base_size = 15, base_rect_size = 1.5) +
          xlab("Group") +
          ylab("log2 Expression") +
          ggtitle(label = "Mean & SEM Plot with Jitter", subtitle = paste(gene.symb, "Probe ID: ", gois[i], sep = "")) +
          stat_summary(fun = "mean", geom = geo, color = "black", size = 0.5, width = 0.2) +
          # Recommended alternative for bar charts in ggplot2:
          #geom_bar(data = sgs, aes(GROUP, MEAN), stat = "identity", color = "black", size = 0.5, width = 0.2) +
          geom_errorbar(data = sgs, aes(GROUP, MEAN, ymin = MEAN - SEM, ymax = MEAN + SEM), size = 1, width = 0.1) + 
          geom_jitter(position = position_jitter(width = 0.1, height = 0, seed = 123), size = 1.5))
      
    } else {
      
      cat("\n")
      stop("Invalid chart.type!\n\n")
    }
    
    printPlots(paste("SingleGene Plot - ", chart.type, " - ", gois[i], sep = ""))
  }
}

#'------------------------------------------------------------------------------
#' @description Quality Control (QC) plots for Agilent array pre-processing.
#'
#' @param exp.mat Expression Matrix whose data are to be plotted.
#' @param targets Agilent "Targets" File.
#' @param stage Pre-processing step label to be used in figure name.
#' @param maplots Number of MA-plots to be drawn.
#'
#' @author //FeAR//
#'------------------------------------------------------------------------------
Agilent.QC.Plots = function(exp.mat, targets, stage = "", maplots = 3)
{
  library(MatrixGenerics) # Enable rowMedians() function
  m = dim(exp.mat)[2]
  
  # Quality Control Boxplots ('las' specifies the style of axis labels)
  graphics::boxplot(exp.mat, names = targets$SampleNumber, las = 1,
          main = paste(stage, "Expression"), xlab = "Sample Index",
          ylab = "log2(intensity)")
  printPlots(paste("Boxplots -", stage))
  
  # Quality Control Density Plots
  plotDensities(exp.mat, main = paste(stage, "Densities"), legend = FALSE)
  printPlots(paste("Densities -", stage))
  
  # Quality Control MA-plots
  # Sampling some (n == maplots) chips
  if (maplots >= m) {
    sample.index = round(seq(from = 1, to = m, length = m))
  } else {
    sample.index = round(seq(from = 1, to = m, length = maplots))
  }
  # 1 chip vs the pseudo-median reference chip
  for (i in sample.index) {
    #oligo::MAplot(exp.mat, which = i, pairs = FALSE, transfo = identity)
    rafalib::maplot(rowMedians(exp.mat[,-i]),
                    exp.mat[,i],
                    xlab = "A (Average log-expression)",
                    ylab = "M (Expression log-ratio)",
                    main = paste0(stage, ": chip ", targets$SampleNumber[i],
                                  " vs pseudo-median reference chip"),
                    n = 5e4,
                    curve.add = TRUE, curve.col = "red3", curve.lwd = 1.5,
                    curve.n = 1e4, pch = 20, cex = 0.1)
    abline(h = c(1,-1), col = "cornflowerblue", lty = 2) # lty = line type
    abline(h = 0, col = "cornflowerblue") # lty = line type
    printPlots(paste0("MA-Plot_", targets$SampleNumber[i], " - ", stage))
  }
}

#'------------------------------------------------------------------------------
#' @description Quality Control (QC) plots for Affymetrix array pre-processing.
#'
#' @param affy.exp Affymetrix object containing expression data.
#' @param exon.flag FALSE for Affy 3'IVT arrays / TRUE  for Affy ST arrays.
#' @param stage Pre-processing step label to be used in figure name.
#' @param maplots Number of MA-plots to be drawn.
#' @param trans.func Function to be used to transform data before boxplotting
#'                   (it can be log2 or identity).
#'
#' @author //FeAR//
#'------------------------------------------------------------------------------
Affymetrix.QC.Plots = function(affy.exp, exon.flag = TRUE, stage = "",
                               maplots = 3, trans.func = identity)
{
  m = dim(affy.exp)["Samples"]
  
  # Quality Control Boxplots
  set.seed(101) # For reproducibility, because oligo::boxplot methods use a
                # sample of the probesets to produce the plot
  main = paste(stage, "Expression")
  if (exon.flag) {
    # For Gene ST and Exon ST arrays, oligo::boxplot() needs the 'target' option,
    # when used with a FeatureSet-like object. 'target' argument can take one of
    # the following values: "probeset", "core" (Gene/Exon), "full" (Exon),
    # "extended" (Exon), describing the summarization target.
    oligo::boxplot(affy.exp, target = "core", names = c(1:m), main = main,
                   xlab = "Sample Index", ylab = "log2(intensity)",
                   transfo = trans.func)
  } else {
    oligo::boxplot(affy.exp, names = c(1:m), main = main, xlab = "Sample Index",
                   ylab = "log2(intensity)", transfo = trans.func)
  }
  printPlots(paste("Boxplots -", stage))
  
  # Quality Control MA-plots
  # Sampling some (n == maplots) chips
  if (maplots >= m) {
    sample.index = round(seq(from = 1, to = m, length = m))
  } else {
    sample.index = round(seq(from = 1, to = m, length = maplots))
  }
  # 1 chip vs the pseudo-median reference chip
  for (i in sample.index) {
    oligo::MAplot(affy.exp, which = i, pairs = FALSE)
    printPlots(paste0("MA-Plot_", i, " - ", stage))
  }
}
