######################################################################################
# sigsQcNegativeControl.R
#
# Implements negative and permutation controls for gene signature QC metrics.
#
# Two functions:
#   .sigsQcNegativeControl() — outer loop dispatcher: iterates over all signatures
#                               in the input list, calling .sigQcNegativeControl for each.
#   .sigQcNegativeControl()  — inner worker: for a single signature, generates random
#                               signatures (negative control) and permuted labels
#                               (permutation control), computes QC metrics via
#                               .compute_without_plots(), and produces summary boxplots.
#
# Negative controls: randomly select the same number of genes from the full expression
#   matrix to form "null" signatures. If QC metrics for the real signature are similar
#   to random signatures, the signature may lack biological coherence.
#
# Permutation controls: keep the same genes but shuffle (permute) gene labels independently
#   for each sample. This destroys inter-gene correlations while preserving marginal
#   distributions, testing whether the observed structure is due to co-expression.
#
# Output structure:
#   out_dir/negative_control/<dataset>/<signature>/sigQC/ — null QC results
#   out_dir/negative_control/<dataset>/<signature>/boxplot_metrics.pdf — comparison plot
#   out_dir/permutation_control/<dataset>/<signature>/sigQC/ — permutation QC results
#   out_dir/permutation_control/<dataset>/<signature>/boxplot_metrics.pdf — comparison plot
#
######################################################################################

######################################################################################
# .sigsQcNegativeControl — Outer loop dispatcher
#
# Iterates over all signatures in genesList and calls .sigQcNegativeControl for each.
#
# @param genesList A list of signatures, each a matrix of gene IDs
# @param expressionMatrixList A list of gene expression matrices (genes x samples)
# @param outputDir Output directory path
# @param studyName Name of the study
# @param numResampling Number of bootstrap resamplings (default 50)
# @param warningsFile Path to warnings log file
# @param logFile Path to log file
#
# @author Andrew Dhawan, Alessandro Barberis
######################################################################################
.sigsQcNegativeControl <- function(genesList, expressionMatrixList, outputDir, studyName, numResampling=50, warningsFile, logFile){
  ###########Check the input
  if(missing(genesList)){
    stop("Need to specify a list of genes. The IDs must match those in the expression matrices.")
  }
  if(missing(expressionMatrixList)){
    #stop("Neet to specify an expression matrix where the rows represent the genes and the columns the samples.")
    stop("Neet to specify a list of expression matrices where the rows represent the genes and the columns the samples.")
  }
  if(missing(outputDir)){
    stop("Need to specify an output directory")
  }
  if(missing(studyName)){
    studyName = "MyStudy"
  }
  if(missing(warningsFile)){
    warningsFile = file.path(outputDir, "warnings.log")
  }
  if(missing(logFile)){
    logFile = file.path(outputDir, "log.log")
  }
  tryCatch({
    re.power = numResampling;
    #Check if the output directory exists. Create otherwise.
    if(!dir.exists(outputDir)){
      dir.create(path=outputDir, showWarnings = F, recursive = T)
    }
    #Open connection to log and warning files
    log.con = file(logFile, open = "a")
    warnings.con = file(warningsFile, open = "a")

    #Loop over the signatures
    for(genes.i in 1:length(genesList)){
      genes = as.matrix(genesList[[genes.i]])
      colnames(genes) = names(genesList)[genes.i]
      #For each signature do a negative control
      # NOTE: Parameter name inconsistency — this function uses numResampling,
      # but .sigQcNegativeControl uses rePower. They refer to the same thing.
      .sigQcNegativeControl(genes, expressionMatrixList, outputDir, studyName, rePower=numResampling, warningsFile, logFile);
    }
  }, error = function(err) {
    cat("", file=log.con, sep="\n")
    cat(paste(Sys.time(),"Errors occurred during in sigQcNegativeControl:", err, sep=" "), file=log.con, sep="\n")
    #stop("Errors occurred during the computation of negative controls")
  }, finally = {
    #cat("---------------------------------------------------------------------------", file=log.con, sep="\n")
    close(log.con)
    close(warnings.con)
  })#END tryCatch
}

######################################################################################
# .sigQcNegativeControl — Per-signature negative and permutation control worker
#
# For a single gene signature:
#   Part 1 (Negative control): Generates rePower random signatures of the same length,
#     computes QC metrics for each, summarizes with quantiles, and plots boxplot comparison.
#   Part 2 (Permutation control): Permutes gene labels within the signature for each sample
#     independently, computes QC metrics on permuted data, and plots boxplot comparison.
#
# @param genes A single-column matrix of gene IDs (the signature)
# @param expressionMatrixList A list of expression matrices
# @param outputDir Output directory path
# @param studyName Study name
# @param rePower Number of resamplings (default 50)
# @param warningsFile Path to warnings log
# @param logFile Path to log
#
# @author Alessandro Barberis
######################################################################################
.sigQcNegativeControl <- function(genes, expressionMatrixList, outputDir, studyName, rePower=50, warningsFile, logFile){
  ###########Check the input
  if(missing(genes)){
    stop("Need to specify a list of genes. The IDs must match those in the expression matrices.")
  }
  if(missing(expressionMatrixList)){
    stop("Neet to specify a list of expression matrices where the rows represent the genes and the columns the samples.")
  }
  if(missing(outputDir)){
    stop("Need to specify an output directory")
  }
  if(missing(studyName)){
    studyName = "MyStudy"
  }
  if(missing(warningsFile)){
    warningsFile = file.path(outputDir, "warnings.log")
  }
  if(missing(logFile)){
    logFile = file.path(outputDir, "log.log")
  }
  tryCatch({
    re.power = rePower;
    #Check if the output directory exists. Create otherwise.
    if(!dir.exists(outputDir)){
      dir.create(path=outputDir, showWarnings = F, recursive = T)
    }
    #Open connection to log and warning files
    log.con = file(logFile, open = "a")
    warnings.con = file(warningsFile, open = "a")

    #Define some useful variables
    len = dim(genes)[1] # Number of genes in the signature

    # ===========================================================================================
    # PART 1: NEGATIVE CONTROLS (random gene selection)
    # ===========================================================================================
    # For each dataset, generate rePower random signatures of the same length as the
    # original signature, compute all QC metrics (via .compute_without_plots), then
    # summarize with quantiles and compare to the original signature's metrics.

    #Loop over datasets
    datasets.num = length(expressionMatrixList)
    datasets.names = names(expressionMatrixList)
    for(dataset.i in 1:datasets.num){
      expressionMatrix = expressionMatrixList[[dataset.i]]
      if(is.null(datasets.names))
        datasetName = paste0("Dataset",dataset.i)
      else
        datasetName = datasets.names[dataset.i]

      data.matrix.ncols = dim(expressionMatrix)[2]
      data.matrix.nrows = dim(expressionMatrix)[1]

      #Define the signature
      gene_sigs_list = list()
      signatureName = colnames(genes)#e.g. "hypoxiaSig"

      # Generate rePower random signatures by randomly selecting `len` genes
      #Loop n times (re-sampling power)
      for(i in 1:re.power){
        # BUG (critical): Uses runif() for discrete gene index selection.
        # runif() returns continuous values in [min, max]. When used as row indices,
        # R truncates to integer, introducing sampling bias:
        #   - Gene at index 1 has near-zero probability (requires runif exactly 1.0)
        #   - Gene at last index gets extra probability from the continuous interval
        #   - Genes can be selected more than once (no uniqueness guarantee)
        # FIX: Replace with sample(1:data.matrix.nrows, size=len, replace=FALSE)
        #Compute a random index vector (fixed: use sample instead of runif for discrete selection)
        random.index.vector = sample(1:data.matrix.nrows, size=len, replace=FALSE);
        #Random signature
        random.genes = as.matrix(rownames(expressionMatrix)[random.index.vector]);
        #Add to the signatures list
        gene_sigs_list[[paste0("NC",i)]] = random.genes;
      }
      names_sigs = names(gene_sigs_list)

      #####sigQC
      cat(paste(Sys.time(), "Computing the Negative Control...", sep=" "), file=log.con, sep="")

      # Run QC pipeline (without plots) on all random signatures at once
      mRNA_expr_matrix = list()
      #names = c(names, gse)
      names = c(datasetName)

      mRNA_expr_matrix[[datasetName]] = expressionMatrix
      sigQC.out_dir = file.path(outputDir, "negative_control", datasetName, signatureName, "sigQC")
      dir.create(path = sigQC.out_dir, showWarnings = FALSE, recursive = T)

      .compute_without_plots(names_datasets = names,
                            gene_sigs_list = gene_sigs_list,
                            names_sigs = names_sigs,
                            mRNA_expr_matrix = mRNA_expr_matrix,
                            out_dir = sigQC.out_dir
                            )


      cat("DONE", file=log.con, sep="\n")

      # Read the aggregated radar chart metrics from all random signatures
      #Read the radarchart_table in the sigQC dir
      metrics.file.path = file.path(sigQC.out_dir, "radarchart_table", "radarchart_table.txt")
      metrics.table = utils::read.table(file = metrics.file.path, header = T, sep = "\t", check.names = F, row.names = 1)

      # Compute summary statistics: mean and quantiles (2.5%, 25%, 50%, 75%, 97.5%)
      # These define the null distribution of each metric
      #Create a matrix with 3 rows (mean, quantiles 0.025 and 0.975)
      neg.controls.summary = matrix(nrow = 6, ncol = dim(metrics.table)[2])
      colnames(neg.controls.summary) = colnames(metrics.table)
      rownames(neg.controls.summary) = c("mean", "Q0.025", "Q0.25", "Q0.5", "Q0.75", "Q0.975")

      #Compute the metrics mean, values for which the 95% of the population fall into
      neg.controls.summary[1,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){mean(t, na.rm=T)})
      neg.controls.summary[2,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.025))})
      neg.controls.summary[3,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.25))})
      neg.controls.summary[4,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.5))})
      neg.controls.summary[5,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.75))})
      neg.controls.summary[6,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.975))})

      # Prepare boxplot overlay: null distribution (gray) + original metric value (red)
      #Define the input variable for the function boxplot.matrix2
      stripchartMatrixList=list()
      stripchartMatrixList[["Negative Control"]]=metrics.table

      # Retrieve the original signature's QC metrics from the main output directory
      #Retrieve data from the original sigQC
      sigQC.out_dir = outputDir
      sig.metrics.file.path = file.path(outputDir, "radarchart_table", "radarchart_table.txt")
      if(file.exists(sig.metrics.file.path)){
        sig.metrics.table = utils::read.table(file = sig.metrics.file.path, header = T, sep = "\t", check.names = F, row.names = 1)
        #TO CHANGE: need to address a table where there could be stored multiple signatures metrics
        #sig.metrics.table = t(sig.metrics.table)
        # Match the row for this specific dataset-signature combination
        # Row names use dots instead of spaces (R's make.names convention)
        #Look for the row label containing the signature name and the current dataset id
       # paste0(gsub(' ','.', names_datasets[i]),'_',gsub(' ','.', names_sigs[k]))
        index = which(rownames(sig.metrics.table) == paste0(gsub(' ','.', datasetName),'_',gsub(' ','.', signatureName))) #which(grepl(paste0(".*(",datasetName,".*",signatureName,"|",signatureName,".*",datasetName,").*"), rownames(sig.metrics.table), ignore.case=TRUE)==T)
        #sig.metrics.table.sigs =

        sig.metrics.table = sig.metrics.table[index,,drop=F]
        stripchartMatrixList[["Original Metric Value"]]=sig.metrics.table
      }

      # Generate boxplot comparing null distribution to original signature metric
      sigQC.out_dir = file.path(outputDir, "negative_control", datasetName, signatureName)
      # Labels for the 14 radar chart metrics
      stripchart_group_names <- c('Relative Med. SD','Skewness',expression(sigma["" >= "10%" ]),expression(sigma["" >= "25%" ]),expression(sigma["" >= "50%" ]),'Coef. of Var.',
                               'Non-NA Prop.','Prop. Expressed',
                               'Autocor.',expression(rho["Mean,Med" ]),
                               expression(rho["PCA1,Med" ]),expression(rho["Mean,PCA1" ]), expression(sigma["PCA1" ]),
                               expression(rho["Med,Z-Med" ]))

      .boxplot.matrix2(x=neg.controls.summary[2:6,], outputDir=sigQC.out_dir, plotName="boxplot_metrics",
                       plotTitle=paste("Boxplot Negative Controls for",signatureName,sep=" "), stripchartMatrixList=stripchartMatrixList,
                       stripchartPch=c(1, 21), stripchartCol = c("gray", "red"), xlab ="Metrics", ylab="Score",group.names=stripchart_group_names)

      #Store the matrix containing the summary of negative controls (i.e. quantiles)
      summary.filePath = file.path(sigQC.out_dir, "neg_controls_summary_table.txt")
      utils::write.table(x = neg.controls.summary,file = summary.filePath,row.names = TRUE, col.names = TRUE,sep=",")

    }#END LOOP OVER DATASETS

    # ===========================================================================================
    # PART 2: PERMUTATION CONTROLS (gene label shuffling)
    # ===========================================================================================
    # For each dataset, keep the same signature genes but permute gene labels independently
    # for each sample. This destroys the inter-gene correlation structure while preserving
    # each gene's marginal expression distribution across samples.
    #
    # The permutation is done per-column (per-sample): for each sample, the expression values
    # of the signature genes are shuffled among themselves. This means different samples get
    # different permutations.
    #
    # ISSUE: The permutation strategy creates an independent random permutation for each
    # sample column. This may be intentional (to fully destroy structure) but it also
    # destroys any sample-level patterns. An alternative approach would be to apply a
    # single permutation to all samples, which would only destroy gene identity while
    # preserving sample-level correlation structure.

   #----the following is for the permutation control----
      #here, we are going to permute the labels of the genes of the signature on all of the samples
      #so looping over each gene signature, we are going to get different dataset matrices such that the rows are permuted
      # we are going to modify the relevant rows of each of the expression matrices for each of the negative controls


    for(dataset.i in 1:datasets.num){
      if(is.null(datasets.names))
        datasetName = paste0("Dataset",dataset.i)
      else
        datasetName = datasets.names[dataset.i]

      expressionMatrix = expressionMatrixList[[dataset.i]]
      signatureName = colnames(genes)#e.g. "hypoxiaSig"

      # Generate rePower permuted expression matrices
      expressionMatrix_perm_list <- list()
      for(i in 1:re.power){
        expressionMatrix_perm <- expressionMatrix
        genes_present <- intersect(rownames(expressionMatrix),genes[,1])
        # Generate a matrix of permutation indices: one independent permutation per sample column
        new_ordering <- replicate(sample(1:length(genes_present)),n=dim(expressionMatrix)[2])
        # Apply per-column permutation to the signature gene rows
        for(col.num in 1:dim(expressionMatrix)[2]){
          expressionMatrix_perm[genes_present,col.num] <- expressionMatrix[genes_present[new_ordering[,col.num]],col.num]
        }
        expressionMatrix_perm_list[[paste0("PC",i)]] <- expressionMatrix_perm
      }

      # Run QC pipeline on all permuted datasets (each permutation becomes a "dataset")
      sigQC.out_dir = file.path(outputDir, "permutation_control", datasetName, signatureName, "sigQC")
      dir.create(path = sigQC.out_dir, showWarnings = FALSE, recursive = T)
      gene_sigs_list <- list()
      gene_sigs_list[[colnames(genes)]] <- genes
      .compute_without_plots(names_datasets = names(expressionMatrix_perm_list),
                            gene_sigs_list = gene_sigs_list,
                            names_sigs = names(gene_sigs_list),
                            mRNA_expr_matrix = expressionMatrix_perm_list,
                            out_dir = sigQC.out_dir
                            )


      cat("DONE", file=log.con, sep="\n")

      # Read and summarize permutation control metrics (same logic as negative control above)
      #Read the radarchart_table in the sigQC dir
      metrics.file.path = file.path(sigQC.out_dir, "radarchart_table", "radarchart_table.txt")
      metrics.table = utils::read.table(file = metrics.file.path, header = T, sep = "\t", check.names = F, row.names = 1)

      #Create a matrix with 3 rows (mean, quantiles 0.025 and 0.975)
      neg.controls.summary = matrix(nrow = 6, ncol = dim(metrics.table)[2])
      colnames(neg.controls.summary) = colnames(metrics.table)
      rownames(neg.controls.summary) = c("mean", "Q0.025", "Q0.25", "Q0.5", "Q0.75", "Q0.975")

      # Compute summary quantiles for the permutation null distribution
      #Compute the metrics mean, values for which the 95% of the population fall into
      neg.controls.summary[1,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){mean(t, na.rm=T)})
      neg.controls.summary[2,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.025))})
      neg.controls.summary[3,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.25))})
      neg.controls.summary[4,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.5))})
      neg.controls.summary[5,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.75))})
      neg.controls.summary[6,] = apply(X=metrics.table,MARGIN = 2,FUN = function(t){stats::quantile(t, na.rm=T, probs=c(0.975))})

      # Prepare boxplot overlay: permutation distribution (gray) + original metric (red)
      #Define the input variable for the function boxplot.matrix2
      stripchartMatrixList=list()
      stripchartMatrixList[["Permutation Control"]]=metrics.table

      #Retrieve data from the original sigQC
      sigQC.out_dir = outputDir
      sig.metrics.file.path = file.path(outputDir, "radarchart_table", "radarchart_table.txt")
      if(file.exists(sig.metrics.file.path)){
        sig.metrics.table = utils::read.table(file = sig.metrics.file.path, header = T, sep = "\t", check.names = F, row.names = 1)
        #TO CHANGE: need to address a table where there could be stored multiple signatures metrics
        #sig.metrics.table = t(sig.metrics.table)
        #Look for the row label containing the signature name and the current dataset id
       # paste0(gsub(' ','.', names_datasets[i]),'_',gsub(' ','.', names_sigs[k]))
        index = which(rownames(sig.metrics.table) == paste0(gsub(' ','.', datasetName),'_',gsub(' ','.', signatureName))) #which(grepl(paste0(".*(",datasetName,".*",signatureName,"|",signatureName,".*",datasetName,").*"), rownames(sig.metrics.table), ignore.case=TRUE)==T)
        #sig.metrics.table.sigs =

        sig.metrics.table = sig.metrics.table[index,,drop=F]
        stripchartMatrixList[["Original Metric Value"]]=sig.metrics.table
      }

      # Generate boxplot for permutation controls
      sigQC.out_dir = file.path(outputDir, "permutation_control", datasetName, signatureName)
      stripchart_group_names <- c('Relative Med. SD','Skewness',expression(sigma["" >= "10%" ]),expression(sigma["" >= "25%" ]),expression(sigma["" >= "50%" ]),'Coef. of Var.',
                               'Non-NA Prop.','Prop. Expressed',
                               'Autocor.',expression(rho["Mean,Med" ]),
                               expression(rho["PCA1,Med" ]),expression(rho["Mean,PCA1" ]), expression(sigma["PCA1" ]),
                               expression(rho["Med,Z-Med" ]))

      .boxplot.matrix2(x=neg.controls.summary[2:6,], outputDir=sigQC.out_dir, plotName="boxplot_metrics",
                       plotTitle=paste("Boxplot Permutation Controls for",signatureName,sep=" "), stripchartMatrixList=stripchartMatrixList,
                       stripchartPch=c(1, 21), stripchartCol = c("gray", "red"), xlab ="Metrics", ylab="Score",group.names=stripchart_group_names)

      #Store the matrix containing the summary of negative controls (i.e. quantiles)
      summary.filePath = file.path(sigQC.out_dir, "perm_controls_summary_table.txt")
      utils::write.table(x = neg.controls.summary,file = summary.filePath,row.names = TRUE, col.names = TRUE,sep=",")

    }

    #----------------------------------------------------

  }, error = function(err) {
    cat("", file=log.con, sep="\n")
    cat(paste(Sys.time(),"Errors occurred during in sigQcNegativeControl:", err, sep=" "), file=log.con, sep="\n")
    #stop("Errors occurred during the computation of negative controls")
  }, finally = {
    #cat("---------------------------------------------------------------------------", file=log.con, sep="\n")
    close(log.con)
    close(warnings.con)
  })#END tryCatch
}
