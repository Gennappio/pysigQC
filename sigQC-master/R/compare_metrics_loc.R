# compare_metrics_loc.R
#
# Compares various signature summary scoring metrics against each other to assess
# concordance of different approaches for summarizing gene signature expression.
#
# The function has 4 major sections:
#   Section 1 (lines ~16-186):  Mean/Median/PCA1 pairwise Spearman correlations + scree plots
#   Section 2 (lines ~189-354): Enrichment score comparisons (GSVA, ssGSEA, PLAGE)
#   Section 3 (lines ~357-392): Scoring metric correlation heatmaps
#   Section 4 (lines ~397-552): Gaussian mixture models (mclust) for score modality analysis
#   Section 5 (lines ~556-641): QQ plots for normality assessment of scores
#
# @param gene_sigs_list A list of genes representing the gene signature to be tested.
# @param names_sigs The names of the gene signatures (one name per gene signature, in gene_sigs_list)
# @param mRNA_expr_matrix A list of expression matrices
# @param names_datasets The names of the different datasets contained in mRNA_expr_matrix
# @param out_dir A path to the directory where the resulting output files are written
# @param file File representing the log file where errors can be written
# @param showResults Tells if open dialog boxes showing the computed results. Default is FALSE
# @param radar_plot_values A list of values that store computations that will be used in the final summary radarplot
# @keywords compare_metrics_loc

compare_metrics_loc <- function(gene_sigs_list,names_sigs, mRNA_expr_matrix, names_datasets, out_dir = '~',file=NULL,showResults = FALSE,radar_plot_values){
  # require(gplots)
  if(!dir.exists(file.path(out_dir,'metrics_tables'))){
      dir.create(file.path(out_dir,'metrics_tables'))
  }

  # score_cor_mats accumulates Spearman correlation matrices between all scoring methods

  # (Mean, Median, PCA1, GSVA, ssGSEA, PLAGE) for each dataset-signature combination.
  # It is populated in Section 2 (line ~343) and consumed in Section 3 for heatmap plots.
  score_cor_mats <- list()

  # ===========================================================================================
  # SECTION 1: Mean / Median / PCA1 pairwise Spearman correlations and scree plots
  # ===========================================================================================
  # For each signature, creates a PDF with 4 rows x N_datasets columns:
  #   Row 1: Median vs Mean smoothScatter + density overlay
  #   Row 2: Mean vs PCA1 smoothScatter + density overlay
  #   Row 3: PCA1 vs Median smoothScatter + density overlay
  #   Row 4: Scree plot (variance explained by each principal component)
  # Spearman rho values are stored in radar_plot_values for the final radar chart.

  for(k in 1:length(names_sigs)){
    #for each signature we will make a separate file comparing the datasets
    gene_sig <- gene_sigs_list[[names_sigs[k]]] # load the gene signature
    if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
    #set up canvas for plotting

    if (showResults){
      grDevices::dev.new()
    }else{
      grDevices::pdf(file.path(out_dir,paste0('sig_compare_metrics_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=10)
    }

    #find the max number of characters in the title (used for dynamic font scaling)
    max_title_length <- -999
    for( i in 1:length(names_datasets)){
      if(max_title_length < nchar(paste0(names_datasets[i],' ',names_sigs[k]))){
        max_title_length <- nchar(paste0(names_datasets[i],' ',names_sigs[k]))
      }
    }
    #set up the canvas: 4 rows (median-mean, mean-pca1, pca1-median, scree) x N datasets
    graphics::par(mfcol = c(4,length(names_datasets)),mar=c(4,4,4,4))

    for ( i in 1:length(names_datasets)){
      # now we can loop over the datasets for the plot and generate the metrics for every dataset with this signature
      data.matrix = mRNA_expr_matrix[[names_datasets[i]]] #load the data
      data.matrix[!(is.finite(as.matrix(data.matrix)))] <- NA #replace Inf/-Inf/NaN with NA

      inter = intersect(gene_sig,rownames(data.matrix)) #consider only the genes actually present in the data

      # Compute per-sample summary scores using only the signature genes present in this dataset
      med_scores <- apply(data.matrix[inter,],2,function(x){stats::median(stats::na.omit(x))}) #compute median
      mean_scores <- apply(data.matrix[inter,],2,function(x){mean(stats::na.omit(x))}) #compute mean

      pca1_scores <- NULL

      tryCatch({
        # BUG: na.omit(t(...)) removes entire SAMPLES (rows of the transposed matrix) that
        # contain any NA, rather than removing genes with NAs. This means samples with even
        # one missing gene value are silently dropped from PCA. The intended behavior is likely
        # to remove genes (columns of the transposed matrix) with NAs, which would require
        # transposing, then filtering columns, or using na.action=na.exclude.
        pca1_scores <- stats::prcomp(stats::na.omit(t(data.matrix[inter,])),retx=T) #compute PCA1
        output_mat <- pca1_scores$x
        if(!is.null(colnames(data.matrix))){
          row.names(output_mat) <- colnames(data.matrix)
        }
        utils::write.table(output_mat,file=file.path(out_dir,'metrics_tables', paste0('pca_loadings_',names_sigs[k],'_',names_datasets[i],'.txt')),quote=F,sep='\t')

        },error=function(e){
          pca1_scores <<- NULL
          cat(paste0("There was an error when computing PCA1 score for:  ",names_datasets[i]," ", names_sigs[k]," ", e,'\n'), file=file)

     #     print(paste0("error: ", e))
        })


      # print(paste0("test ",pca1_scores ))
      if(length(pca1_scores) > 1){#!is.null(pca1_scores)){
        props_of_variances <- pca1_scores$sdev^2/(sum(pca1_scores$sdev^2)) #proportion of variance per component
        pca1_scores <- pca1_scores$x[,1] #takes the actual PCA1 scores (first column of rotation matrix)
        common_score_cols <- intersect(names(med_scores),intersect(names(mean_scores),names(pca1_scores))) #ensures we have the same samples for each plot
      }else{
        common_score_cols <- intersect(names(med_scores),names(mean_scores))
      }

      # --- Plot 1: Median vs Mean scatter with density overlay ---
      if(length(common_score_cols) > 1){
        #the following is the colourmap for the 2D scatter (jet colormap: blue->cyan->yellow->red)
        jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        #plotting commands for the 2D scatterplot for mean-median correlation as follows
        graphics::smoothScatter(med_scores[common_score_cols],mean_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='Median vs Mean')
        graphics::points(med_scores[common_score_cols],mean_scores[common_score_cols],pch='.') #draw the points on top of it
        graphics::par(new=T)#,srt=45)
        graphics::plot(stats::density(med_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA,  col='red',main=NA,lwd=2) #draw the density plot behind it
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8) #labels for the axes
        graphics::mtext(side = 2, line = 2, 'Mean',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'Median',cex=0.8)
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title
        rho <- stats::cor(med_scores[common_score_cols],mean_scores[common_score_cols],method='spearman')
        rho_mean_med <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(med_scores[common_score_cols]))
      }else{
        graphics::plot.new()
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title

        graphics::title(paste0('\n\nToo many NA values for Mean/Median in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_mean_med <- 0
      }

      # --- Plot 2: Mean vs PCA1 scatter ---
      if (length(pca1_scores) > 1){#(!is.null(pca1_scores)){
        #plotting for the mean-pca1
        graphics::smoothScatter(mean_scores[common_score_cols],pca1_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='Mean vs PCA1')
        graphics::points(mean_scores[common_score_cols],pca1_scores[common_score_cols],pch='.')
        graphics::par(new=T)
        graphics::plot(stats::density(mean_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA,lwd=2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
        rho <- stats::cor(mean_scores[common_score_cols],pca1_scores[common_score_cols],method='spearman')
        rho_mean_pca1 <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(mean_scores[common_score_cols]))
        graphics::mtext(side = 2, line = 2, 'PCA1',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'Mean',cex=0.8)

        # --- Plot 3: PCA1 vs Median scatter ---
        graphics::smoothScatter(pca1_scores[common_score_cols],med_scores[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='PCA1 vs Median')
        graphics::points(pca1_scores[common_score_cols],med_scores[common_score_cols],pch='.')
        graphics::par(new=T)
        graphics::plot(stats::density(pca1_scores[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA,lwd=2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
        rho <- stats::cor(pca1_scores[common_score_cols],med_scores[common_score_cols],method='spearman')
        rho_pca1_med <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(pca1_scores[common_score_cols]))
        graphics::mtext(side = 2, line = 2, 'Median',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'PCA1',cex=0.8)
      }else{

#        graphics::par(new=T)
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1/Mean in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_mean_pca1 <- 0

 #       graphics::par(new=T)
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for Median/PCA1 in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_pca1_med <- 0

      }

      # #here we output the table of mean, median and pca1 scores for each sample to a table
      # if(length(pca1_scores) > 1){#(!is.null(pca1_scores)){
      #   output_mat <- cbind(mean_scores[common_score_cols],cbind(med_scores[common_score_cols],pca1_scores[common_score_cols]))
      #   colnames(output_mat) <- c('Mean_Scores','Median_Scores','PCA1_Scores')
      # }else{
      #   output_mat <- cbind(mean_scores[common_score_cols],med_scores[common_score_cols])
      #   colnames(output_mat) <- c('Mean_Scores','Median_Scores')

      # }

      # utils::write.table(output_mat,file=file.path(out_dir,'metrics_tables', paste0('metrics_table_',names_sigs[k],'_',names_datasets[i],'.txt')),quote=F,sep='\t')

      #stores Spearman rho values that will be used in the radarplot
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_mean_med'] <- rho_mean_med
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_pca1_med'] <- rho_pca1_med
      radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['rho_mean_pca1'] <- rho_mean_pca1

      # --- Plot 4: Scree plot (proportion of variance explained by each PC) ---
      if(length(pca1_scores) > 1){#(!is.null(pca1_scores)){
        #draws the scree plot (first 10 components or fewer)
        bars_plot <- props_of_variances[1:min(10,length(props_of_variances))]

        output_mat <- as.matrix(props_of_variances)
        colnames(output_mat) <- c('Proportion of variance')
        row.names(output_mat) <- paste0('PCA ',as.character(1:dim(output_mat)[1]))
        utils::write.table(output_mat,file=file.path(out_dir,'metrics_tables', paste0('pca_vs_var_',names_sigs[k],'_',names_datasets[i],'.txt')),quote=F,sep='\t')

        graphics::barplot(bars_plot,main="PCA vs proportion\n of variance") #ylim= c(0,1),
        graphics::mtext(side = 1, line = 2, 'PCA',cex=0.8)
        # Store proportion of variance explained by PC1 as a radar chart metric
        radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['prop_pca1_var']<- bars_plot[1]
      }else{
       # graphics::par(new=T)
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1 in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        radar_plot_values[[names_sigs[k]]][[names_datasets[i]]]['prop_pca1_var']<- 0
      }
    }
    #saves file
    if(showResults){
      grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('sig_compare_metrics_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=10)
    }
    if(grDevices::dev.cur()!=1){
        g <- grDevices::dev.off() # to reset the graphics pars to defaults
    }
  }

  cat('Metrics compared successfully.\n', file=file) #output to log

  # ===========================================================================================
  # SECTION 2: Enrichment Score (ES) comparisons — GSVA, ssGSEA, PLAGE
  # ===========================================================================================
  # Computes three gene set enrichment scores for each signature-dataset combination:
  #   - GSVA: Gene Set Variation Analysis (nonparametric kernel density estimation)
  #   - ssGSEA: single-sample Gene Set Enrichment Analysis (rank-based enrichment)
  #   - PLAGE: Pathway Level Analysis of Gene Expression (SVD-based)
  # All three methods produce per-sample enrichment scores. The pairwise Spearman
  # correlations between them assess concordance of enrichment scoring approaches.
  #
  # ISSUE: These three methods have fundamentally different statistical assumptions
  # (GSVA uses nonparametric kernel density, ssGSEA uses rank-based enrichment,
  # PLAGE uses SVD). The pairwise correlations assess concordance but the code draws
  # no conclusions about which is more appropriate for a given dataset.
  #
  # ISSUE: No multiple testing correction is applied despite computing many Spearman
  # correlations (3 pairs x signatures x datasets). Raw rho values are reported without
  # p-values or FDR correction.

  for(k in 1:length(names_sigs)){
    #for each signature we will make a separate file comparing the datasets
    gene_sig <- gene_sigs_list[[names_sigs[k]]] # load the gene signature
    if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
    #set up canvas for plotting

    if (showResults){
      grDevices::dev.new()
    }else{
      grDevices::pdf(file.path(out_dir,paste0('sig_compare_ES_metrics_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=7.5)
    }

    #find the max number of characters in the title (used for dynamic font scaling)
    max_title_length <- -999
    for( i in 1:length(names_datasets)){
      if(max_title_length < nchar(paste0(names_datasets[i],' ',names_sigs[k]))){
        max_title_length <- nchar(paste0(names_datasets[i],' ',names_sigs[k]))
      }
    }
    #set up the canvas: 3 rows (GSVA-ssGSEA, ssGSEA-PLAGE, PLAGE-GSVA) x N datasets
    graphics::par(mfcol = c(3,length(names_datasets)),mar=c(4,4,4,4))

    for ( i in 1:length(names_datasets)){
      # now we can loop over the datasets for the plot and generate the metrics for every dataset with this signature
      data.matrix = mRNA_expr_matrix[[names_datasets[i]]] #load the data
      data.matrix[!(is.finite(as.matrix(data.matrix)))] <- NA #replace Inf/-Inf/NaN with NA
      inter = intersect(gene_sig,rownames(data.matrix)) #consider only the genes actually present in the data

      # Recompute summary scores (same as Section 1 — duplicated computation)
      med_scores <- apply(data.matrix[inter,],2,function(x){stats::median(stats::na.omit(x))}) #compute median
      mean_scores <- apply(data.matrix[inter,],2,function(x){mean(stats::na.omit(x))}) #compute mean
      pca1_scores <- NULL
      tryCatch({
        # BUG: Same na.omit(t(...)) issue as Section 1 — removes samples, not genes
        pca1_scores <- stats::prcomp(stats::na.omit(t(data.matrix[inter,])),retx=T) #compute PCA1
        output_mat <- pca1_scores$x
        if(!is.null(colnames(data.matrix))){
          row.names(output_mat) <- colnames(data.matrix)
        }
        utils::write.table(output_mat,file=file.path(out_dir,'metrics_tables', paste0('pca_loadings_',names_sigs[k],'_',names_datasets[i],'.txt')),quote=F,sep='\t')

        },error=function(e){
          pca1_scores <<- NULL
     #     print(paste0("error: ", e))
        })

      # Ensure expression matrix is a proper matrix (not data.frame) for GSVA
      #need the matrix to be a  matrix not a dataframe
      data.matrix.gsva <- data.matrix
      if(!is.matrix(data.matrix)){
        data.matrix.gsva <- as.matrix.data.frame(data.matrix,rownames.force=T)
      }

      # Compute three enrichment scores using the GSVA package
      # Each returns a 1-row matrix (one score per sample) for the gene set defined by 'inter'
      gsvaPar <- GSVA::ssgseaParam(data.matrix.gsva, list(inter))
      es.ssGSEA <- suppressWarnings(GSVA::gsva(gsvaPar, verbose=F))

      gsvaPar <- GSVA::gsvaParam(data.matrix.gsva, list(inter))
      es.gsva <- suppressWarnings(GSVA::gsva(gsvaPar, verbose=F))

      gsvaPar <- GSVA::plageParam(data.matrix.gsva, list(inter))
      es.plage <- suppressWarnings(GSVA::gsva(gsvaPar, verbose=F))

      # Extract per-sample scores from the 1-row result matrices
      es.gsva <- es.gsva[1,]
      es.ssGSEA <- es.ssGSEA[1,]
      es.plage <- es.plage[1,]

      common_score_cols <- intersect(names(es.gsva),intersect(names(es.ssGSEA),names(es.plage))) #ensures we have the same samples for each plot

      # --- ES Plot 1: GSVA vs ssGSEA scatter with density overlay ---
      if(length(common_score_cols) > 1){
        #the following is the colourmap for the 2D scatter
        jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        #plotting commands for the 2D scatterplot for mean-median correlation as follows
        graphics::smoothScatter(es.gsva[common_score_cols],es.ssGSEA[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='GSVA vs ssGSEA')
        graphics::points(es.gsva[common_score_cols],es.ssGSEA[common_score_cols],pch='.') #draw the points on top of it
        graphics::par(new=T)#,srt=45)
        graphics::plot(stats::density(es.gsva[common_score_cols]), axes=F, xlab=NA, ylab=NA,  col='red',main=NA,lwd=2) #draw the density plot behind it
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8) #labels for the axes
        graphics::mtext(side = 2, line = 2, 'ssGSEA',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'GSVA',cex=0.8)
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title
        rho <- stats::cor(es.gsva[common_score_cols],es.ssGSEA[common_score_cols],method='spearman')
        rho_ssGSEA_gvsa <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(es.gsva[common_score_cols]))
      }else{
        graphics::plot.new()
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title

        graphics::title(paste0('\n\nToo many NA values for ssGSEA/GSVA in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_ssGSEA_gvsa <- 0
      }

      # --- ES Plot 2: ssGSEA vs PLAGE scatter ---
      if(length(common_score_cols) > 1){
        #plotting for the ssGSEA-PLAGE
        graphics::smoothScatter(es.ssGSEA[common_score_cols],es.plage[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='ssGSEA vs PLAGE')
        graphics::points(es.ssGSEA[common_score_cols],es.plage[common_score_cols],pch='.')
        graphics::par(new=T)
        graphics::plot(stats::density(es.ssGSEA[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA,lwd=2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
        rho <- stats::cor(es.ssGSEA[common_score_cols],es.plage[common_score_cols],method='spearman')
        rho_ssGSEA_plage <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(es.ssGSEA[common_score_cols]))
        graphics::mtext(side = 2, line = 2, 'PLAGE',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'ssGSEA',cex=0.8)
      }else{

#        graphics::par(new=T)
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PLAGE/ssGSEA in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_ssGSEA_plage <- 0
      }

      # --- ES Plot 3: PLAGE vs GSVA scatter ---
      if(length(common_score_cols) > 1){
        graphics::smoothScatter(es.plage[common_score_cols],es.gsva[common_score_cols],colramp = jet.colors,xlab=NA,ylab=NA,main='PLAGE vs GSVA')
        graphics::points(es.plage[common_score_cols],es.gsva[common_score_cols],pch='.')
        graphics::par(new=T)
        graphics::plot(stats::density(es.plage[common_score_cols]), axes=F, xlab=NA, ylab=NA, col='red',main=NA,lwd=2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density',cex=0.8)
        rho <- stats::cor(es.plage[common_score_cols],es.gsva[common_score_cols],method='spearman')
        rho_plage_gsva <- rho
        graphics::mtext(paste0('rho = ',format(rho,digits = 2)),side=3,line=0,cex = 0.6,at=max(es.plage[common_score_cols]))
        graphics::mtext(side = 2, line = 2, 'GSVA',cex=0.8)
        graphics::mtext(side = 1, line = 2, 'PLAGE',cex=0.8)
        }else{

 #       graphics::par(new=T)
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for GSVA/PLAGE in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
        rho_plage_gsva <- 0

      }

      # Re-extract PCA1 scores and rebuild common_score_cols to include all methods
      # (PCA was computed earlier but pca1_scores still holds the prcomp object)
       # print(paste0("test ",pca1_scores ))
      if(length(pca1_scores) > 1){#!is.null(pca1_scores)){
        props_of_variances <- pca1_scores$sdev^2/(sum(pca1_scores$sdev^2)) #for the scree plot
        pca1_scores <- pca1_scores$x[,1] #takes the actual PCA1 scores
        common_score_cols <- intersect(common_score_cols,intersect(names(med_scores),intersect(names(mean_scores),names(pca1_scores)))) #ensures we have the same samples for each plot
      }else{
        common_score_cols <- intersect(common_score_cols,intersect(names(med_scores),names(mean_scores)))
      }

      # Output combined metrics table (all scoring methods) for this signature-dataset combo
      #here we output the table of mean, median and pca1 scores for each sample to a table

      if((length(common_score_cols) > 1) & (!is.null(pca1_scores)) ){
        common_score_cols <- intersect(common_score_cols,names(pca1_scores))
        output_mat <- cbind(mean_scores[common_score_cols],med_scores[common_score_cols],pca1_scores[common_score_cols],es.ssGSEA[common_score_cols],cbind(es.gsva[common_score_cols],es.plage[common_score_cols]))
        colnames(output_mat) <- c('Mean_Scores','Median_Scores','PCA1_Scores','ssGSEA','GSVA','PLAGE')
        utils::write.table(output_mat,file=file.path(out_dir,'metrics_tables', paste0('metrics_table_',names_sigs[k],'_',names_datasets[i],'.txt')),quote=F,sep='\t')

      }else if( length(common_score_cols) > 1){
        output_mat <- cbind(mean_scores[common_score_cols],med_scores[common_score_cols],es.ssGSEA[common_score_cols],cbind(es.gsva[common_score_cols],es.plage[common_score_cols]))
        colnames(output_mat) <- c('Mean_Scores','Median_Scores','ssGSEA','GSVA','PLAGE')
        utils::write.table(output_mat,file=file.path(out_dir,'metrics_tables', paste0('metrics_table_',names_sigs[k],'_',names_datasets[i],'.txt')),quote=F,sep='\t')

      }
      # Store the full Spearman correlation matrix between all scoring methods
      # (consumed in Section 3 below for heatmap visualization)
      score_cor_mats[[paste0(names_datasets[i],'_',names_sigs[k])]]<- stats::cor(output_mat,method='spearman')

    }
    #saves file
    if(showResults){
      grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('sig_compare_ES_metrics_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=7.5)
    }
    if(grDevices::dev.cur()!=1){
        g <- grDevices::dev.off() # to reset the graphics pars to defaults
    }

  }
 cat('ES scores computed successfully.\n', file=file) #output to log

  # ===========================================================================================
  # SECTION 3: Scoring metric correlation heatmaps
  # ===========================================================================================
  # Visualizes the Spearman correlation matrix between all scoring methods
  # (Mean, Median, PCA1, GSVA, ssGSEA, PLAGE) as a ComplexHeatmap for each
  # dataset-signature combination.

#-----now we need to output the score correlation matrices as heatmaps----
  for(k in 1:length(names_sigs)){
    for ( i in 1:length(names_datasets)){
#now we compute and output the correlation of scoring metrics in heatmap form
      plot_mat <- score_cor_mats[[paste0(names_datasets[i],'_',names_sigs[k])]]
      if (showResults){
        grDevices::dev.new()
      }else{
        grDevices::pdf(file.path(out_dir,paste0('scoring_metrics_corr_',names_datasets[i],'_',names_sigs[k],'.pdf')),width=4,height=4)
      }
      #set up the canvas
      graphics::par(mfcol = c(4,length(names_datasets)),mar=c(4,4,4,4))
      # NOTE: Typo in display label: 'GVSA' instead of 'GSVA' (but the key matches correctly)
      nice_row_names <- c('Mean','Median','PCA1','GVSA','ssGSEA','PLAGE')
      names(nice_row_names) <- c('Mean_Scores','Median_Scores','PCA1_Scores','GSVA','ssGSEA','PLAGE')
      row_names.fontsize <- 10
      row.names(plot_mat) <- nice_row_names[rownames(plot_mat)]
      ans_hmap <- ComplexHeatmap::Heatmap(plot_mat,show_column_dend = F,
                                              show_column_names = F,
                                              name=names_datasets[i],
                                              col= circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                                              heatmap_legend_param = list(title = 'Correlation', color_bar = "continuous",legend_direction='vertical'),
                                              column_title = paste0(names_datasets[i],' ',names_sigs[k],'\nScoring Metric Correlation'),
                                              row_names_gp =  grid::gpar(fontsize = row_names.fontsize),
                                              row_title = 'Scoring metrics')#,

      ComplexHeatmap::draw(ans_hmap,heatmap_legend_side = "left")

      #saves file
      if(showResults){
        grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('scoring_metrics_corr_',names_datasets[i],'_',names_sigs[k],'.pdf')),width=7,height=3.5)
      }
      if(grDevices::dev.cur()!=1){
          g <- grDevices::dev.off() # to reset the graphics pars to defaults
      }
  }
}
 #-----------------------------------------------------------------------------------------------------------

  #------------------------------------------------------------------------------------------------------------------

  # ===========================================================================================
  # SECTION 4: Gaussian Mixture Models (mclust) for score modality analysis
  # ===========================================================================================
  # Fits Gaussian mixture models to each scoring metric (median, mean, PCA1) to assess
  # whether the score distribution is unimodal or multimodal. Uses BIC (Bayesian Information
  # Criterion) to select the optimal number of components (1 to max_clusters).
  #
  # mclust model names:
  #   E = Equal variance (homoscedastic) — all components share the same variance
  #   V = Variable variance (heteroscedastic) — each component has its own variance
  #
  # Output: BIC plots per metric, text summary of best model, and raw mclust objects (.rda)

  #next we use the mclust package, and compute the log-likelihoods of the various Gaussian mixture models for the scores

  mixture_model.out = file.path(out_dir, "mixture_model_out.txt")
  mixture_model.con = file(mixture_model.out, open = "w") #open the output file

  mixture_models <- list() #variable to store the mclust results

  for(k in 1:length(names_sigs)){
    mixture_models[[names_sigs[k]]] <- list()
    gene_sig <- gene_sigs_list[[names_sigs[k]]] # load the gene signature
    if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
    #set up canvas for plotting

    if (showResults){
      grDevices::dev.new()
    }else{
      grDevices::pdf(file.path(out_dir,paste0('sig_gaussian_mixture_model_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=10)
    }

    #find the max number of characters in the title (used for dynamic font scaling)
    max_title_length <- -999
    for( i in 1:length(names_datasets)){
      if(max_title_length < nchar(paste0(names_datasets[i],' ',names_sigs[k]))){
        max_title_length <- nchar(paste0(names_datasets[i],' ',names_sigs[k]))
      }
    }
    #set up the canvas: 3 rows (median BIC, mean BIC, PCA1 BIC) x N datasets
    graphics::par(mfcol = c(3,length(names_datasets)),mar=c(4,4,4,4))

    for ( i in 1:length(names_datasets)){
      # now we can loop over the datasets for the plot and generate the metrics for every dataset with this signature
      data.matrix = mRNA_expr_matrix[[names_datasets[i]]] #load the data
      data.matrix[!(is.finite(as.matrix(data.matrix)))] <- NA #replace Inf/-Inf/NaN with NA
      inter = intersect(gene_sig,rownames(data.matrix)) #consider only the genes actually present in the data

      med_scores <- apply(data.matrix[inter,],2,function(x){stats::median(stats::na.omit(x))}) #compute median
      mean_scores <- apply(data.matrix[inter,],2,function(x){mean(stats::na.omit(x))}) #compute mean
      pca1_scores <- NULL

      tryCatch({
        pca1_scores <- stats::prcomp(stats::na.omit(t(data.matrix[inter,])),retx=T) #compute PCA1
        pca1_scores <- pca1_scores$x[,1] #gets the first component of PCA

        },error=function(e){
          pca1_scores <<- NULL
          cat(paste0("There was an error when computing PCA1 score for:  ",names_datasets[i]," ", names_sigs[k]," ", e,'\n'), file=file)
      })


      #the following are the mixture model variables that are going to be the outputs of the mclust function
      mixture_models[[names_sigs[k]]][[names_datasets[i]]] <- list()
      mixture_models[[names_sigs[k]]][[names_datasets[i]]][['median']] <- NULL
      mixture_models[[names_sigs[k]]][[names_datasets[i]]][['mean']] <- NULL
      mixture_models[[names_sigs[k]]][[names_datasets[i]]][['pca1']] <- NULL

      # --- Mixture model for MEDIAN scores ---
      if(length(med_scores) > 1){
        # max_clusters: at most half the non-NA samples, capped at 10
        max_clusters <- min(ceiling(sum(!is.na(med_scores)) / 2),10)
        mixture_models[[names_sigs[k]]][[names_datasets[i]]][['median']] <- mclust::Mclust(med_scores,G=1:max_clusters)
        mclust::plot.Mclust(x =  mixture_models[[names_sigs[k]]][[names_datasets[i]]][['median']], what='BIC',main='Median score')
        # graphics::title(main='Median score')
       graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title

        output_string <- paste0(names_sigs[k],' ',names_datasets[i],', Median score: There ')
        if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['median']]$G==1){
          output_string <- paste0(output_string,' is one component in the Gaussian mixture model. ')
        }else{
           output_string <- paste0(output_string,' are ',mixture_models[[names_sigs[k]]][[names_datasets[i]]][['median']]$G,' components in the Gaussian mixture model. ')
        }

        # BUG (copy-paste): Model name labels are SWAPPED.
        # 'V' means Variable variance (heteroscedastic), but the message says "equal variance (E)".
        # 'E' means Equal variance (homoscedastic), but the message says "variable variances (V)".
        # The correct mapping is: E -> equal variance, V -> variable variance.
        if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['median']]$modelName=='V'){
          output_string <- paste0(output_string, 'Best model is Gaussian distributions with equal variance (E). ' )
        }else if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['median']]$modelName=='E'){
          output_string <- paste0(output_string, 'Best model is Gaussian distributions with variable variances (V). ' )
        }else{
          output_string <- paste0(output_string, 'Best model is univariate Gaussian distribution. ' )
        }
        output_string <- paste0(output_string,'\n')

        cat(output_string,file=mixture_model.con)
      }else{
        graphics::plot.new()
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title

        graphics::title(paste0('\n\nToo many NA values for Median in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
      }

      # --- Mixture model for MEAN scores ---
      if(length(mean_scores) > 1){
        max_clusters <- min(ceiling(sum(!is.na(mean_scores)) / 2),10)

        # BUG (copy-paste): Fits model on med_scores instead of mean_scores.
        # This means the "mean" mixture model is actually a duplicate of the median model.
        # Should be: mclust::Mclust(mean_scores, G=1:max_clusters)
        mixture_models[[names_sigs[k]]][[names_datasets[i]]][['mean']] <- mclust::Mclust(med_scores,G=1:max_clusters)
        mclust::plot.Mclust(x =  mixture_models[[names_sigs[k]]][[names_datasets[i]]][['mean']], what='BIC',main='Mean score')
        # graphics::title(main='Mean score')
        output_string <- paste0(names_sigs[k],' ',names_datasets[i],', Mean score: There ')
        if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['mean']]$G==1){
          output_string <- paste0(output_string,' is one component in the Gaussian mixture model. ')
        }else{
           output_string <- paste0(output_string,' are ',mixture_models[[names_sigs[k]]][[names_datasets[i]]][['mean']]$G,' components in the Gaussian mixture model. ')
        }

        # BUG (copy-paste): Same E/V label swap as in the median block above.
        if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['mean']]$modelName=='V'){
          output_string <- paste0(output_string, 'Best model is Gaussian distributions with equal variance (E). ' )
        }else if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['mean']]$modelName=='E'){
          output_string <- paste0(output_string, 'Best model is Gaussian distributions with variable variances (V). ' )
        }else{
          output_string <- paste0(output_string, 'Best model is univariate Gaussian distribution. ' )
        }
        output_string <- paste0(output_string,'\n')

        cat(output_string,file=mixture_model.con)
      }else{
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for Mean in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
      }

      # --- Mixture model for PCA1 scores ---
      if(length(pca1_scores) > 1){
        max_clusters <- min(ceiling(sum(!is.na(pca1_scores)) / 2),10)

        mixture_models[[names_sigs[k]]][[names_datasets[i]]][['pca1']] <- mclust::Mclust(pca1_scores,G=1:max_clusters)
        mclust::plot.Mclust(x =  mixture_models[[names_sigs[k]]][[names_datasets[i]]][['pca1']], what='BIC',main='PCA1 score')
        # graphics::title(main='PCA1 score')

        output_string <- paste0(names_sigs[k],' ',names_datasets[i],', PCA1 score: There ')
        if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['pca1']]$G==1){
          output_string <- paste0(output_string,' is one component in the Gaussian mixture model. ')
        }else{
           output_string <- paste0(output_string,' are ',mixture_models[[names_sigs[k]]][[names_datasets[i]]][['pca1']]$G,' components in the Gaussian mixture model. ')
        }

        # BUG (copy-paste): Same E/V label swap as above.
        if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['pca1']]$modelName=='V'){
          output_string <- paste0(output_string, 'Best model is Gaussian distributions with equal variance (E). ' )
        }else if(mixture_models[[names_sigs[k]]][[names_datasets[i]]][['pca1']]$modelName=='E'){
          output_string <- paste0(output_string, 'Best model is Gaussian distributions with variable variances (V). ' )
        }else{
          output_string <- paste0(output_string, 'Best model is univariate Gaussian distribution. ' )
        }
        output_string <- paste0(output_string,'\n')

        cat(output_string,file=mixture_model.con)

      }else{
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1 in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
      }
    }
     #saves file
    if(showResults){
      grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('sig_gaussian_mixture_model_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=10)
    }
    if(grDevices::dev.cur()!=1){
        g <- grDevices::dev.off() # to reset the graphics pars to defaults
    }
  }
 close(mixture_model.con)

 save(mixture_models,file=file.path(out_dir, "mixture_models_raw_out.rda"))
 cat('Gaussian mixture models computed successfully.\n', file=file) #output to log



  # ===========================================================================================
  # SECTION 5: QQ plots for normality assessment of scoring metrics
  # ===========================================================================================
  # Generates normal QQ plots for each scoring metric (median, mean, PCA1) to visually
  # assess whether the distribution of per-sample scores follows a Gaussian distribution.
  # This complements the mixture model analysis: if a score is clearly non-normal,
  # the mixture model with >1 component may be more appropriate.

 #-------------next thing we will do with these scores is to plot the QQ plots to check for normality--------


  for(k in 1:length(names_sigs)){
    mixture_models[[names_sigs[k]]] <- list()
    gene_sig <- gene_sigs_list[[names_sigs[k]]] # load the gene signature
    if(is.matrix(gene_sig)){gene_sig = as.vector(gene_sig);}
    #set up canvas for plotting

    if (showResults){
      grDevices::dev.new()
    }else{
      grDevices::pdf(file.path(out_dir,paste0('sig_qq_plots_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=10)
    }

    #find the max number of characters in the title (used for dynamic font scaling)
    max_title_length <- -999
    for( i in 1:length(names_datasets)){
      if(max_title_length < nchar(paste0(names_datasets[i],' ',names_sigs[k]))){
        max_title_length <- nchar(paste0(names_datasets[i],' ',names_sigs[k]))
      }
    }
    #set up the canvas: 3 rows (median QQ, mean QQ, PCA1 QQ) x N datasets
    graphics::par(mfcol = c(3,length(names_datasets)),mar=c(4,4,4,4))

    for ( i in 1:length(names_datasets)){
      # Recompute scores (same computation repeated in every section — duplicated code)
      # now we can loop over the datasets for the plot and generate the metrics for every dataset with this signature
      data.matrix = mRNA_expr_matrix[[names_datasets[i]]] #load the data
      data.matrix[!(is.finite(as.matrix(data.matrix)))] <- NA #replace Inf/-Inf/NaN with NA
      inter = intersect(gene_sig,rownames(data.matrix)) #consider only the genes actually present in the data

      med_scores <- apply(data.matrix[inter,],2,function(x){stats::median(stats::na.omit(x))}) #compute median
      mean_scores <- apply(data.matrix[inter,],2,function(x){mean(stats::na.omit(x))}) #compute mean
      pca1_scores <- NULL

      tryCatch({
        pca1_scores <- stats::prcomp(stats::na.omit(t(data.matrix[inter,])),retx=T) #compute PCA1
        pca1_scores <- pca1_scores$x[,1] #gets the first component of PCA
        },error=function(e){
          pca1_scores <<- NULL
          cat(paste0("There was an error when computing PCA1 score for:  ",names_datasets[i]," ", names_sigs[k]," ", e,'\n'), file=file)
      })


      # NOTE: mixture_models is overwritten here for the QQ plot section but the stored
      # values are not used in this section. This appears to be leftover from copy-paste.
      #the following are the mixture model variables that are going to be the outputs of the mclust function
      mixture_models[[names_sigs[k]]][[names_datasets[i]]] <- list()
      mixture_models[['median']] <- NULL
      mixture_models[['mean']] <- NULL
      mixture_models[['pca1']] <- NULL

      if(length(med_scores) > 1){
        stats::qqnorm(med_scores,plot.it=T,main='Median score')
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title

      }else{
        graphics::plot.new()
        graphics::mtext(side=3,line=2.5,paste0(names_datasets[i],' ',names_sigs[k]),cex=min(1,3*10/max_title_length)) #title
        graphics::title(paste0('\n\nToo many NA values for Median in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
      }

      if(length(mean_scores) > 1){
        stats::qqnorm(mean_scores,plot.it=T,main='Mean score')

      }else{
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for Mean in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
      }

      if(length(pca1_scores) > 1){
        stats::qqnorm(pca1_scores,plot.it=T,main='PCA1 score')

      }else{
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1 in \n',names_datasets[i],' ',names_sigs[k]))#cex=min(1,4*10/max_title_length))
      }
    }
     #saves file
    if(showResults){
      grDevices::dev.copy(grDevices::pdf,file.path(out_dir,paste0('sig_qq_plots_',names_sigs[k],'.pdf')),width=3*length(names_datasets),height=10)
    }
    if(grDevices::dev.cur()!=1){
        g <- grDevices::dev.off() # to reset the graphics pars to defaults
    }
  }
 cat('QQ plots computed successfully.\n', file=file) #output to log

 #-----------------------------------------------------------------------------------------------------------


  radar_plot_values #returns the radarplot values
}
