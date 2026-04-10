# eval_struct_loc.R
#
# Evaluates gene signature structure via hierarchical clustering heatmaps and
# biclustering analysis. Produces:
#   Section 1: Expression heatmaps with hierarchical clustering (ComplexHeatmap)
#              with optional sample annotation overlays (covariates)
#   Section 2: Biclustering on z-transformed + binarized data using BCCC algorithm,
#              displayed as heatmaps with continuous color scale
#   Section 3: Same biclustering but displayed on the binarized (0/1) heatmap
#
# The biclustering sections (2 and 3) contain heavily duplicated code — the same
# z-transform, binarization, and BCCC computation is repeated 4 times total
# (twice for checking if any biclusters exist, twice for plotting).
#
# @param gene_sigs_list A list of genes representing the gene signature to be tested.
# @param names_sigs The names of the gene signatures (one name per gene signature, in gene_sigs_list)
# @param mRNA_expr_matrix A list of expression matrices
# @param names_datasets The names of the different datasets contained in mRNA_expr_matrix
# @param covariates A list containing a sub-list of 'annotations' and 'colors' which contains the annotation matrix for the given dataset and the associated colours with which to plot in the expression heatmap
# @param out_dir A path to the directory where the resulting output files are written
# @param file File representing the log file where errors can be written
# @param showResults Tells if open dialog boxes showing the computed results. Default is FALSE
# @param radar_plot_values A list of values that store computations that will be used in the final summary radarplot
# @keywords eval_struct_loc

eval_struct_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix,names_datasets,covariates, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  # library(gplots)
  # require(biclust)
  # require(ComplexHeatmap)
  # par(cex.main=0.7,cex.lab = 0.6,oma=c(2,2,3,2),mar=c(2,2,2,2))

  #find the max number of characters in the title for fontsize purposes
  max_title_length <- -999
  for(k in 1:length(names_sigs)){
    for( i in 1:length(names_datasets)){
      if(max_title_length < nchar(paste0(names_datasets[i],' ',names_sigs[k]))){
        max_title_length <- nchar(paste0(names_datasets[i],' ',names_sigs[k]))
      }
    }
  }

  # ===========================================================================================
  # SECTION 1: Hierarchical clustering heatmaps
  # ===========================================================================================
  # For each signature, creates one PDF per dataset showing the expression heatmap
  # of signature genes with hierarchical clustering on samples (columns).
  # If covariates are provided, they are displayed as annotation bars above the heatmap.

  # --- Pre-processing: ensure all heatmaps for a given signature have the same gene rows ---
  # Genes missing from a dataset are padded with min(expression) values.
  #
  # BUG (statistical): Padding missing genes with the dataset minimum expression value
  # artificially introduces perfectly correlated rows (identical values across all samples).
  # This distorts both hierarchical clustering distances and any downstream autocorrelation.
  # A better approach would be to pad with NA and handle missing values in clustering.

  #first we do a check to ensure that all heatmpas have the same genes as rows and add NA values if not
  #first let's get all possible rows expressed of the signature across all samples
  all_row_names <- list() #stores the union of gene names present across all datasets for each signature
  for (k in 1:length(names_sigs)){
    all_row_names[[names_sigs[k]]] <- c()
    gene_sig <- gene_sigs_list[[names_sigs[k]]] #load the signature
    if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
    for (i in 1:length(names_datasets)){

      data.matrix = mRNA_expr_matrix[[names_datasets[i]]] #load the dataset
      inter <- intersect(gene_sig, row.names(data.matrix)) #consider only the genes present in the dataset

      all_row_names[[names_sigs[k]]] <- c(all_row_names[[names_sigs[k]]],inter)
    }
     all_row_names[[names_sigs[k]]] <- unique(all_row_names[[names_sigs[k]]])
  }

  # Build padded expression matrices so all datasets have the same gene rows per signature
  #now that we know the unique rownames for each signature, we should go through and generate the heatmap matrices with appended min values
  sig_scores_all_mats <- list() #this will be the list of signature score matrices
  for (k in 1:length(names_sigs)){
    sig_scores_all_mats[[names_sigs[k]]] <- list()
    gene_sig <- gene_sigs_list[[names_sigs[k]]] #load the signature
    if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
    for (i in 1:length(names_datasets)){
      data.matrix = mRNA_expr_matrix[[names_datasets[i]]] #load the dataset
      inter <- intersect(gene_sig, row.names(data.matrix)) #consider only the genes present in the dataset
      sig_scores <- (as.matrix(data.matrix[inter,])) #extract the signature expression submatrix
      #now let's see how many rows of NA values we need to add
      rows_needed <- setdiff(all_row_names[[names_sigs[k]]],inter)
      # BUG (operator precedence): `length(rows_needed >0)` evaluates `rows_needed > 0` first
      # (comparing character vector to numeric 0, which coerces and always returns a logical vector
      # of the same length), then takes length() of that logical vector.
      # Should be: `length(rows_needed) > 0`
      # The condition happens to work correctly in practice because length(logical_vector) > 0
      # whenever rows_needed is non-empty, but it is true for the wrong reason.
      # FIXED:
      if(length(rows_needed) > 0){
        # BUG (statistical): Padding with min(sig_scores) creates rows with identical values,
        # which have zero variance and produce artificial correlation patterns.
        # FIX: Pad with NA instead and use na.rm=TRUE downstream
        sig_scores <- rbind(sig_scores,matrix(NA,nrow=length(rows_needed),ncol=dim(sig_scores)[2]))
        row.names(sig_scores) <- c(inter,rows_needed)
      }
      sig_scores_all_mats[[names_sigs[k]]][[names_datasets[i]]] <- sig_scores#[gene_sig[,1],]
    }
  }

  # --- Generate heatmaps using ComplexHeatmap ---
  for(k in 1:length(names_sigs)){
    gene_sig <- gene_sigs_list[[names_sigs[k]]] #load the signature
    if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
    # Build a list of ComplexHeatmap objects (one per dataset) using sapply
    hmaps <- sapply(1:length(names_datasets),function(i) {
      #this is a subroutine to create a list of heatmaps stored in hmaps that can be plotted

      sig_scores <- sig_scores_all_mats[[names_sigs[k]]][[names_datasets[i]]]
      tryCatch({
        if (length(covariates[[names_datasets[i]]]) ==0){
          # --- Case 1: No covariates — plain heatmap ---

          #here we will set the parameters for the font size etc in the heatmaps
          dim.pdf = dim(sig_scores); #size of the heatmap matrix
          h = dim.pdf[1]; #number of genes (rows)
          if(h<20){ #if less than 20 genes, use standard font; otherwise scale down
            row_names.fontsize = 12
          }else{
            row_names.fontsize=5/log10(h) # logarithmic font scaling for many genes
          }
          #the following creates the heatmap (column dendrogram disabled, row dendrogram enabled by default)
          ans_hmap <- ComplexHeatmap::Heatmap((sig_scores),show_column_dend = F,
                                              show_column_names = F,
                                              name=names_datasets[i],
                                              heatmap_legend_param = list(title = names_datasets[i], color_bar = "continuous",legend_direction='vertical'),
                                              column_title = paste0(names_datasets[i]),
                                              row_names_gp =  grid::gpar(fontsize = row_names.fontsize),
                                              row_title = 'Genes')#,
        }else{
          # --- Case 2: With covariates — heatmap with annotation bars ---
          if (is.vector(covariates[[names_datasets[i]]][['annotations']])){
            # Annotations provided as a named vector (one annotation per sample)
            #this is the case if the user has provided further annotations that they may want to plot alongside the top of the heatmap
            ha1 = ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(covariates[[names_datasets[i]]][['annotations']][intersect(names(covariates[[names_datasets[i]]][['annotations']]),colnames(sig_scores))]),
                                                    col=covariates[[names_datasets[i]]][['colors']],
                                                    na_col="grey")#,
            #show_annotation_name = TRUE)#, col = list(type = c("a" = "red", "b" = "blue")
          }else{
            # Annotations provided as a matrix/data.frame (multiple annotations per sample)
            ha1 = ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(covariates[[names_datasets[i]]][['annotations']][intersect(rownames(covariates[[names_datasets[i]]][['annotations']]),colnames(sig_scores)),]),
                                                    col=covariates[[names_datasets[i]]][['colors']],
                                                    na_col="grey",
                                                    which="column")#,
            #show_annotation_name = TRUE)#, col = list(type = c("a" = "red", "b" = "blue"),
          }
          #this is for setting fontsize of rownames of heatmap
          dim.pdf = dim(sig_scores);
          h = dim.pdf[1];
          if(h<20){
            row_names.fontsize = 12
          }else{
            row_names.fontsize=5/log10(h)
          }

          # row_names_gp = grid::gpar(fontsize = row_names.fontsize)
          #create the heatmap with annotation bars on top
          ans_hmap <- ComplexHeatmap::Heatmap((sig_scores),show_column_dend = F,
                                              show_column_names = F,
                                              name=names_datasets[i],
                                              heatmap_legend_param = list(title = names_datasets[i], color_bar = "continuous",legend_direction='vertical'),
                                              column_title = paste0(names_datasets[i]),
                                              row_names_gp = grid::gpar(fontsize = row_names.fontsize),
                                              row_title = 'Genes', top_annotation=ha1)#,

        }
      },
      error=function(err){
        #print(paste0("There was an error, likely due to NA values in ",names[i] ," : ", err))
        cat(paste0('Error when creating expression heatmaps for ',names_datasets[i],' ', names_sigs[k],': ',err,'\n'), file=file) #output to llog file

        graphics::plot.new()
        graphics::title(paste0('\n \n \n',names_datasets[i],' ',names_sigs[k])) #makes a graphics object
      })
      ans_hmap #return the heatmap to be added to the list for plotting

      #grab_grob()
    })

    # Concatenate heatmaps into a HeatmapList and draw one PDF per dataset,
    # using that dataset's clustering to order all heatmaps.
    #the following loop concatenates the list of heatmaps to be plotted into one heatmaplist (complexheatmap object), and then
    #cycles through each of the datasets for the clustering, plotting each case of heatmap
    for ( i in 1:length(names_datasets)){
      all_hmaps <- hmaps[[1]]

      if (length(hmaps) > 1){
        tmp <- lapply(hmaps[2:length(hmaps)],function(x) all_hmaps <<- ComplexHeatmap::add_heatmap(all_hmaps,x))
      }
      #creates graphic device
      if (showResults){
        grDevices::dev.new()
      } else{
        grDevices::pdf(file.path(out_dir, paste0('sig_eval_struct_clustering_',names_datasets[i],'_',names_sigs[k],'.pdf')),width=5*(length(names_datasets)),height=10)#paste0(out_dir,'/sig_autocor_hmps.pdf'))
      }
      # Draw the combined heatmap list, using the i-th dataset's dendrogram for column ordering
      #the following plots the heatmap list with the clustering done on the current dataset
      ComplexHeatmap::draw(all_hmaps,heatmap_legend_side = "left",annotation_legend_side = "left", main_heatmap = names_datasets[i])
      #saves the pdf
      if(showResults){
        grDevices::dev.copy(grDevices::pdf,file.path(out_dir, paste0('sig_eval_struct_clustering_',names_datasets[i],'_',names_sigs[k],'.pdf')),width=5*(length(names_datasets)),height=10)#paste0(out_dir,'/sig_autocor_hmps.pdf'))
      }
      cat('Expression heatmaps saved successfully.\n', file=file) #ouptuts to log file
      if(grDevices::dev.cur()!=1){
          g<- grDevices::dev.off() #closes graphics device
        }
      }
  }
  # if(grDevices::dev.cur()!=1){
  #   g<- grDevices::dev.off() #closes graphics device
  # }

  # draw.heatmaps(hmaps,names)

  # ===========================================================================================
  # SECTION 2: Biclustering on z-transformed data (continuous heatmap display)
  # ===========================================================================================
  # Uses the BCCC (Bicluster Constant Column Correlation) algorithm from the biclust package
  # to find biclusters — subsets of genes and samples that show correlated expression patterns.
  #
  # Algorithm: Cheng & Church (2000) biclustering
  #   delta: maximum allowed mean squared residue score (higher = looser clusters)
  #   alpha: scaling factor for row/column deletion (higher = more aggressive pruning)
  #   number: maximum number of biclusters to find
  #
  # ISSUE: BCCC parameters (delta=1, alpha=1.5, number=50) appear to be arbitrary defaults
  # with no domain-specific justification. The original BCXmotifs call (commented out) used
  # different parameters. No sensitivity analysis is provided.
  #
  # ISSUE: num_cols_chosen is computed but NEVER USED — it was a parameter for the
  # commented-out BCXmotifs method. This is dead code.
  #
  # ISSUE (binarization threshold): After z-transformation, the threshold is computed as
  # min + (max-min)/2 = midpoint of the range. For z-scored data, a more principled threshold
  # would be 0 (the mean of z-scores) or use a data-driven method like Otsu's thresholding.
  # The midpoint depends heavily on outliers and has no statistical meaning.
  #
  # BUG (division by zero in z-transform): If a gene has zero variance (constant expression
  # across samples), sd() returns 0 and the division produces NaN/Inf. These propagate
  # through na.omit() and silently change the effective number of genes.

  # --- First pass: check if ANY biclusters exist across all dataset/signature combos ---
  #---------before we do the biclustering plotting, let's determine whether there are ANY biclusters among the datasets/sig combo

  save_bicluster <- FALSE

  for (i in 1:(length(names_datasets)* length(names_sigs))){
      #here we compute all the biclusters

      # Convert linear index i to 2D grid indices (dataset_ind, sig_ind)
      #first convert the index i to an array index for the grid size length(names_datasets) by length(names_sigs)
      dataset_ind <- i %% length(names_datasets)
      if (dataset_ind == 0 ){
        dataset_ind <- length(names_datasets)
      }
      sig_ind <- ceiling(i/length(names_datasets))
      gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]

      if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}

      data.matrix = mRNA_expr_matrix[[names_datasets[dataset_ind]]] #load the data

      inter = intersect(gene_sig, row.names(data.matrix)) #only consider the genes present in the dataset

      sig_scores <- as.matrix(data.matrix[inter,])
      sig_scores[!is.finite(sig_scores)] <- NA #make sure the scores aren't infinte

      # Z-transform each gene independently (center and scale across samples)
      #need to standardize here the matrix before binarizing
      for (gene in inter){
        gene_values <- as.numeric(sig_scores[gene,])
        gene_sd <- stats::sd(gene_values, na.rm=TRUE)
        if (is.na(gene_sd) || gene_sd == 0) {
          sig_scores[gene,] <- 0
        } else {
          sig_scores[gene,] <- (gene_values - mean(gene_values, na.rm=TRUE)) / gene_sd
        }
      }

      # Binarize: values <= midpoint threshold become 0, values > threshold become 1
      #the following does the binarization of the matrix
      threshold <- min(stats::na.omit(t(sig_scores)))+(max(stats::na.omit(t(sig_scores)))-min(stats::na.omit(t(sig_scores))))/2
      x <- stats::na.omit(t(sig_scores)) #> threshold) * 1
      x[x<=threshold] <- 0
      x[x>threshold] <- 1

      # ISSUE (dead code): num_cols_chosen was used by the commented-out BCXmotifs method.
      # The current BCCC call does not use this parameter at all.
      #the following if statement creates the parameters for the biclustering algorithm; namely the number of columns to take in each
      #repetition of the biclustering algorithm (we use BCCC) - this is dependent on the size of the matrix

      if (dim(sig_scores)[2] > 40) {
        num_cols_chosen <- 20
      }else if (dim(sig_scores)[2] > 20){
        num_cols_chosen <- 10
      }else if (dim(sig_scores)[2] > 10){
        num_cols_chosen <- 5
      }else{
        num_cols_chosen <- 2
      }

      # Xmotif <- biclust(x, method=BCXmotifs(), number=50, alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))
      Xmotif <- biclust::biclust(x, method=biclust::BCCC(), delta=1,alpha=1.5, number=50)# alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))
      #the above performs the biclustering

      #if more than 1 bicluster then save_bicluster is true, otherwise false
      if(Xmotif@Number > 1){
        save_bicluster = TRUE
      }
    }

  #-----------------------------------------------------------------------------------------------------

  # --- Second pass: if any biclusters found, regenerate and plot them with continuous heatmap ---
  # NOTE: This is a near-exact duplicate of the first pass above. The biclustering is
  # recomputed from scratch rather than cached from the first pass.

  if(save_bicluster){
    #creates new plot
    grDevices::dev.new()

    graphics::par(cex.main=0.8,cex.lab = 0.8,oma=c(4,2,2,2),mar=c(4,4,4,4)) #sets graphics parameters

    hmaps <- lapply(1:(length(names_datasets)* length(names_sigs)),function(i) {
      #here we define a list of heatmaps for binarized biclustering data
      # Convert linear index to 2D grid indices
      #first convert the index i to an array index for the grid size length(names_datasets) by length(names_sigs)
      dataset_ind <- i %% length(names_datasets)
      if (dataset_ind == 0 ){
        dataset_ind <- length(names_datasets)
      }
      sig_ind <- ceiling(i/length(names_datasets))
      gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]
      if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
      data.matrix = mRNA_expr_matrix[[names_datasets[dataset_ind]]] #load the data
      inter = intersect(gene_sig, row.names(data.matrix)) #only consider the genes present in the dataset
      sig_scores <- as.matrix(data.matrix[inter,])
      sig_scores[!is.finite(sig_scores)] <- NA #make sure the scores aren't infinte

      # Z-transform (duplicated from first pass)
      #need to standardize here the matrix before binarizing
      for (gene in inter){
        gene_values <- as.numeric(sig_scores[gene,])
        gene_sd <- stats::sd(gene_values, na.rm=TRUE)
        if (is.na(gene_sd) || gene_sd == 0) {
          sig_scores[gene,] <- 0
        } else {
          sig_scores[gene,] <- (gene_values - mean(gene_values, na.rm=TRUE)) / gene_sd
        }
      }

      # Binarize (duplicated from first pass)
      #the following does the binarization of the matrix
       threshold <- min(stats::na.omit(t(sig_scores)))+(max(stats::na.omit(t(sig_scores)))-min(stats::na.omit(t(sig_scores))))/2
      x <- stats::na.omit(t(sig_scores)) #> threshold) * 1
      x[x<=threshold] <- 0
      x[x>threshold] <- 1
      # print(paste0('thresh ',threshold))
      # x <- biclust::binarize(stats::na.omit(t(sig_scores)))#discretize(stats::na.omit(t(sig_scores)),nof=10,quant=F)

      # ISSUE (dead code): num_cols_chosen not used by BCCC (same as above)
      #the following if statement creates the parameters for the biclustering algorithm; namely the number of columns to take in each
      #repetition of the biclustering algorithm (we use BCCC) - this is dependent on the size of the matrix

      if (dim(sig_scores)[2] > 40) {
        num_cols_chosen <- 20
      }else if (dim(sig_scores)[2] > 20){
        num_cols_chosen <- 10

      }else if (dim(sig_scores)[2] > 10){
        num_cols_chosen <- 5
      }else{
        num_cols_chosen <- 2
      }

      # Xmotif <- biclust(x, method=BCXmotifs(), number=50, alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))
      Xmotif <- biclust::biclust(x, method=biclust::BCCC(), delta=1,alpha=1.5, number=50)# alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))
      #the above performs the biclustering
      # Plot: if >1 bicluster, show heatmap with CONTINUOUS (z-transformed) values + bicluster overlay
      #if more than 1 bicluster then we output the heatmap with the bicluster on it
      #otherwise we just output the empty title
      if(Xmotif@Number > 1){
        biclust::heatmapBC(stats::na.omit(t(sig_scores)),bicResult=Xmotif,col = gplots::colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
        #	biclust::heatmapBC(x,bicResult=Xmotif,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')

        graphics::title(paste0('\n \n \nBivariate clustering\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]),cex=min(1,4*10/max_title_length))
        graphics::axis(1,at=1:length(rownames(sig_scores)), labels=rownames(sig_scores),las=2,tck=0,cex.axis=0.6)
        #mtext(rownames(sig_scores))
        # }else if (Xmotif@Number == 1){
        # 	biclust::heatmapBC(stats::na.omit(t(sig_scores)),bicResult=Xmotif,number=1,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
      }else{
        #print(paste0("Zero or one co-clusters found: ", names[i]))
        graphics::plot.new()
        graphics::title(paste0('\n\n\n <=1 bivariate clusters for\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]),cex=min(1,4*10/max_title_length))
        cat(paste0('<= 1 bi-clusters found for: ', names_datasets[dataset_ind],' ',names_sigs[sig_ind],'\n'), file=file)

      }
      grab_grob() #grab the heatmap or empty title and add it to the list of heatmaps that will be plotted
    })

    draw.heatmaps(hmaps,names_datasets,names_sigs) #this draws the heatmaps list in a grid
    #saves the biclustering datas
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_eval_bivariate_clustering.pdf'),width=4*(length(names_datasets)),height=4*(length(names_sigs)))
    cat('Bi-clustering completed successfully\n', file=file)

    if(grDevices::dev.cur()!=1){
      g<- grDevices::dev.off() #closes graphics device
    }
  }else{
    cat(paste0('Bi-clustering completed successfully. No bi-clusters found among signature/dataset combinations.\n'),file=file)
  }

  # ===========================================================================================
  # SECTION 3: Biclustering on binarized data (binary 0/1 heatmap display)
  # ===========================================================================================
  # Same biclustering as Section 2, but the heatmap displays the BINARIZED matrix
  # instead of the continuous z-transformed values. This shows the exact 0/1 pattern
  # that was input to the BCCC algorithm.
  #
  # NOTE: This entire section is a near-exact duplicate of Section 2. The only difference
  # is that heatmapBC is called with `x` (binarized) instead of `na.omit(t(sig_scores))`
  # (continuous z-transformed).

 #----------------------------------------------------------------------------------------------
 #---------before we do the binarized biclustering plotting, let's determine whether there are ANY biclusters among the datasets/sig combo

 # --- First pass (again): check for biclusters (duplicated from Section 2) ---
 save_bicluster = FALSE
 for (i in 1:(length(names_datasets) * length(names_sigs))){

    #convert index i into an array index for the grid of heatmaps
    dataset_ind <- i %% length(names_datasets)
    if (dataset_ind == 0 ){
      dataset_ind <- length(names_datasets)
    }
    sig_ind <- ceiling(i/length(names_datasets))
    gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]
    if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
    data.matrix = mRNA_expr_matrix[[names_datasets[dataset_ind]]] #load the data
    inter = intersect(gene_sig, row.names(data.matrix)) #consider only genes present in both cases

    sig_scores <- as.matrix(data.matrix[inter,])
    sig_scores[!is.finite(sig_scores)] <- NA #make sure the data is finite
    # Z-transform (duplicated)
    #standardise by z-transform
    for (gene in inter){
      gene_values <- as.numeric(sig_scores[gene,])
      gene_sd <- stats::sd(gene_values, na.rm=TRUE)
      if (is.na(gene_sd) || gene_sd == 0) {
        sig_scores[gene,] <- 0
      } else {
        sig_scores[gene,] <- (gene_values - mean(gene_values, na.rm=TRUE)) / gene_sd
      }
    }
    # Binarize (duplicated)
    #binarize the matrix
    threshold <- min(stats::na.omit(t(sig_scores)))+(max(stats::na.omit(t(sig_scores)))-min(stats::na.omit(t(sig_scores))))/2
    x <- stats::na.omit(t(sig_scores)) #> threshold) * 1
    x[x<=threshold] <- 0
    x[x>threshold] <- 1
    # x <- biclust::binarize(stats::na.omit(t(sig_scores)))#discretize(stats::na.omit(t(sig_scores)),nof=10,quant=F)
    # ISSUE (dead code): num_cols_chosen not used by BCCC
    #the following if statement helps to decide the parameters for the BCCC algorithm for biclustering

    if (dim(sig_scores)[2] > 40) {
      num_cols_chosen <- 20
    }else if (dim(sig_scores)[2] > 20){
      num_cols_chosen <- 10
    }else if (dim(sig_scores)[2] > 10){
      num_cols_chosen <- 5
    }else{
      num_cols_chosen <- 2
    }
    # following performs the biclustering
    # Xmotif <- biclust(x, method=BCXmotifs(), number=50, alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))
    Xmotif <- biclust::biclust(x, method=biclust::BCCC(), delta=1,alpha=1.5, number=50)
    #if there is more then 1 bicluster then will output the heatmap, otherwise not
    if(Xmotif@Number > 1){
      save_bicluster = TRUE

    }
  }
 #----------------------------------------------------------------------------------------------

  # --- Second pass: plot biclusters with BINARIZED heatmap display ---
  #same thing, but with binarized heatmaps (not the exact heatmaps) underneath the clustering
  #create a new graphics object
  if(save_bicluster){
    grDevices::dev.new()

    # pdf(file.path(out_dir,'sig_eval_bivariate_clustering_binarized_maps.pdf'),width=10,height=10)

    graphics::par(cex.main=0.8,cex.lab = 0.8,oma=c(4,2,2,2),mar=c(4,4,4,4)) #set graphics parameters
    hmaps <- lapply(1:(length(names_datasets) * length(names_sigs)),function(i) {

      #convert index i into an array index for the grid of heatmaps
      dataset_ind <- i %% length(names_datasets)
      if (dataset_ind == 0 ){
        dataset_ind <- length(names_datasets)
      }
      sig_ind <- ceiling(i/length(names_datasets))
      gene_sig <- gene_sigs_list[[names_sigs[sig_ind]]]
      if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
      data.matrix = mRNA_expr_matrix[[names_datasets[dataset_ind]]] #load the data
      inter = intersect(gene_sig, row.names(data.matrix)) #consider only genes present in both cases

      sig_scores <- as.matrix(data.matrix[inter,])
      sig_scores[!is.finite(sig_scores)] <- NA #make sure the data is finite
      # Z-transform (duplicated)
      #standardise by z-transform
      for (gene in inter){
        gene_values <- as.numeric(sig_scores[gene,])
        gene_sd <- stats::sd(gene_values, na.rm=TRUE)
        if (is.na(gene_sd) || gene_sd == 0) {
          sig_scores[gene,] <- 0
        } else {
          sig_scores[gene,] <- (gene_values - mean(gene_values, na.rm=TRUE)) / gene_sd
        }
      }
      # Binarize (duplicated)
      #binarize the matrix
      threshold <- min(stats::na.omit(t(sig_scores)))+(max(stats::na.omit(t(sig_scores)))-min(stats::na.omit(t(sig_scores))))/2
      x <- stats::na.omit(t(sig_scores)) #> threshold) * 1
      x[x<=threshold] <- 0
      x[x>threshold] <- 1
      # x <- biclust::binarize(stats::na.omit(t(sig_scores)))#discretize(stats::na.omit(t(sig_scores)),nof=10,quant=F)
      # ISSUE (dead code): num_cols_chosen not used by BCCC
      #the following if statement helps to decide the parameters for the BCCC algorithm for biclustering

      if (dim(sig_scores)[2] > 40) {
        num_cols_chosen <- 20
      }else if (dim(sig_scores)[2] > 20){
        num_cols_chosen <- 10
      }else if (dim(sig_scores)[2] > 10){
        num_cols_chosen <- 5
      }else{
        num_cols_chosen <- 2
      }
      # following performs the biclustering
      # Xmotif <- biclust(x, method=BCXmotifs(), number=50, alpha=0.5, nd=floor(dim(sig_scores)/num_cols_chosen), ns=num_cols_chosen, sd=floor(dim(sig_scores)/num_cols_chosen))
      Xmotif <- biclust::biclust(x, method=biclust::BCCC(), delta=1,alpha=1.5, number=50)
      # Plot: if >1 bicluster, show heatmap with BINARIZED (0/1) values + bicluster overlay
      #if there is more then 1 bicluster then will output the heatmap, otherwise will output an empty plot with title only
      if(Xmotif@Number > 1){
        #biclust::heatmapBC(stats::na.omit(t(sig_scores)),bicResult=Xmotif,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
        biclust::heatmapBC(x,bicResult=Xmotif,col = gplots::colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')

        graphics::title(paste0('\n \n \nBivariate clustering\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]),cex=min(1,4*10/max_title_length))
        graphics::axis(1,at=1:length(rownames(sig_scores)), labels=rownames(sig_scores),las=2,tck=0,cex.axis=0.6)
        #mtext(rownames(sig_scores))
        # }else if (Xmotif@Number == 1){
        # 	biclust::heatmapBC(stats::na.omit(t(sig_scores)),bicResult=Xmotif,number=1,col = colorpanel(100,"blue","white","red"), xlab='Gene ID',ylab='Sample')
      }else{
        # print(paste0("Zero or one co-clusters found: ", names[i]))
        graphics::plot.new()
        graphics::title(paste0('\n\n\n <=1 bivariate clusters for\n',names_datasets[dataset_ind],' ',names_sigs[sig_ind]),cex=min(1,4*10/max_title_length))
      }
      grab_grob() #outputs the graphics object to the list
    })

    draw.heatmaps(hmaps,names_datasets,names_sigs) #draws the heatmaps
    # the following saves the file
    grDevices::dev.copy(grDevices::pdf,file.path(out_dir,'sig_eval_bivariate_clustering_binarized_maps.pdf'),width=4*(length(names_datasets)),height=4*(length(names_sigs)))
    if(grDevices::dev.cur()!=1){
      g<- grDevices::dev.off() #closes graphics device
    }
  }else{
    cat(paste0('Binarized bi-clustering completed successfully. No binarized bi-clusters found among signature/dataset combinations.\n'),file=file)
  }
  if(grDevices::dev.cur()!=1){
      g<- grDevices::dev.off() #closes graphics device
  }
  radar_plot_values #returns the values for the final plotting function
}
