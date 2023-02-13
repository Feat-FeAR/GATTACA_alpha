# Header Info ------------------------------------------------------------------
# STALKER Functions
#
# Special functions used in GATTACA_alpha
#
# a //FeAR// R-script - 13-Feb-2023
#



#'------------------------------------------------------------------------------
#' @description Save a graphical output to `<folderPrefix> Figures`
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
#' @description Create the annotation data frame. NOTE: this function has been
#'              evolved to `array_create_annot()` function in `cmatools` package.
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
#'              label. NOTE: this function has been evolved to `descriptives()`
#'              function in `cmatools` package.
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
#' @description Plot single gene comparison chart. NOTE: this function has been
#'              evolved to `quick_chart()` function in `cmatools` package.
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
      gene.symb = paste(ann[gois[i], grepl("Symbol", colnames(ann))],
                        " - ", sep = "")
    }
    
    if (chart.type == "BP") {
      
      # NOTICE: When in a for loop, you have to explicitly print your resulting
      # ggplot object
      print(
        ggplot(data = sgex, aes(var.groups, var.expr)) +
          theme_bw(base_size = 15, base_rect_size = 1.5) +
          xlab("Group") +
          ylab("log2 Expression") +
          ggtitle(label = "Box Plot with Jitter",
                  subtitle = paste(gene.symb,"Probe ID: ", gois[i], sep = "")) +
          geom_boxplot(width = 0.5, size = 0.5, notch = TRUE,
                       outlier.shape = NA) +
          stat_summary(fun = "mean", geom = geo, color = "red3", size = 2) +
          geom_jitter(position = position_jitter(width = 0.1, height = 0,
                                                 seed = 123), size = 1.5))
      
    } else if (chart.type == "BC" | chart.type == "MS") {
      
      print(
        ggplot(data = sgex, aes(var.groups, var.expr)) +
          theme_bw(base_size = 15, base_rect_size = 1.5) +
          xlab("Group") +
          ylab("log2 Expression") +
          ggtitle(label = "Mean & SEM Plot with Jitter",
                  subtitle = paste(gene.symb, "Probe ID: ", gois[i], sep = "")) +
          stat_summary(fun = "mean", geom = geo, color = "black", size = 0.5,
                       width = 0.2) +
          geom_errorbar(data = sgs, aes(GROUP, MEAN, ymin = MEAN - SEM,
                                        ymax = MEAN + SEM),
                        size = 1, width = 0.1) + 
          geom_jitter(position = position_jitter(width = 0.1, height = 0,
                                                 seed = 123), size = 1.5))
      
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
