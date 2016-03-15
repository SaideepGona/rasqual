#' Calculate sample-specic offsets for RASQUAL
#'
#'  Calculate offset values for each gene in each condition that
#'  correct for library size and GC content bias.  
#'
#' @param counts Matrix of read counts (genes in rows, samples in columns)
#' @param gene_metadata Matrix of gene metadata (required columns: gene_id, precentage_gc_content). 
#' Used for GC content bias correction. 
#' @param gc_correct If true then correct for GC content bias in addition to library size.
#'
#' @return Matrix of gene-specifc offsets.
#' @export
rasqualCalculateSampleOffsets <- function(counts, gene_metadata = NULL, method = "library_size", gc_correct = TRUE){
  
  if(gc_correct == TRUE){
    #Make sure that gene_metadata has the required columns
    assertthat::assert_that(assertthat::has_name(gene_metadata, "gene_id"))
    assertthat::assert_that(assertthat::has_name(gene_metadata, "percentage_gc_content"))
  }
  
  #Calculate library sizes
  library_size = colSums(counts)
  size_factors = library_size/mean(library_size) #Standardise
  size_matrix = matrix(rep(size_factors, nrow(counts)), nrow = nrow(counts), byrow = TRUE)
  rownames(size_matrix) = rownames(counts)
  colnames(size_matrix) = colnames(counts)
  
  #Apply GC correction
  if(gc_correct == TRUE){
    gc_factor = rasqualGcCorrection(counts, gene_metadata)
    size_matrix = size_matrix * gc_factor
  }
  return(size_matrix)
}

#' Helper function for rasqualGcCorrection
#' 
#' @author Natsuhiko Kumasaka
Quantile <- function(x,k=20){
  x=rank(x,ties="random")
  z=rep(0,length(x))
  for(i in 1:k){
    z = z+as.numeric(x<=quantile(x,i/k,na.rm=T))
  }
  k-z
}

#' Helper function for rasqualGcCorrection
#' 
#' @author Natsuhiko Kumasaka
randomize <- function(x,g=NULL){
    if(is.null(g)){
      n=ncol(x);
      t(apply(x,1,function(xx){xx[order(runif(n))]}))
    }else{
      for(i in unique(g)){
        x[,g==i]=randomize(x[,g==i,drop=F])
      }
      x
    }
  }

#' Estimate the effect of GC-bias for each feature in each sample.
#' 
#' This function does not correct for differences in library size between samples.
#'
#' @param Y Matrix of read counts
#' @param gene_metadata Data frame with gene metadata.
#' (required columnss: gene_id, percentage_gc_content)
#' @param PLOT 
#'
#' @return Matrix of GC-bias offsets for each gene in each condition.
#' @export
#' @author Natsuhiko Kumasaka
rasqualGcCorrection <- function(Y,gene_metadata,PLOT=F){
  
  #Extract GC vector from the metadata matrix
  gene_metadata = as.data.frame(gene_metadata)
  rownames(gene_metadata) = gene_metadata$gene_id
  gene_metadata = gene_metadata[rownames(Y),]
  gcvec = gene_metadata$percentage_gc_content
  
  #Perfrorm GC correction
  bin=Quantile(gcvec,200);
  x=sort(unlist(lapply(split(gcvec,bin),mean)))
  S=apply(Y,2,function(y){unlist(lapply(split(y,bin),sum))[as.character(0:199)]});
  Fs=log(t(t(S)/apply(S,2,sum))/apply(S,1,sum)*sum(as.numeric(S)));
  Gs=apply(Fs,2,function(y){smooth.spline(x,y,spar=1)$y}); 
  if(PLOT){
    par(mfcol=c(5,5),mar=c(2,2,2,2)); 
    for(i in 1:ncol(Y)){
      plot(Fs[,i])
      lines(Gs[,i],col=2)
    }
    matplot(x,Gs,type="l",col=2,lty=1)
  }
  result_matrix = exp(Gs[bin+1,])
  #Add gene ids to rows
  rownames(result_matrix) = rownames(Y)
  return(result_matrix)
}

rasqualMakeCovariates <- function(counts, size_factors) {
  
  #Map parameters to Natsuhiko's variables
  Y = counts
  K = size_factors
  n=ncol(Y)
  
  # fpm calculation
  fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9
  
  # Singular value decomposition
  fpkm.svd   = svd((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd))
  fpkm.svd.r = svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))
  
  # Covariate selection
  sf=log(apply(Y,2,sum))
  covs=fpkm.svd$v[,1:sum(fpkm.svd$d[-n]>fpkm.svd.r$d[-n])]
  if(cor(sf,covs[,1])^2<0.9){covs=cbind(sf, covs)}
  
  # Write covariates
  return(covs)
}

