#' Convert a data frame into a GRanges object
#' 
#' Seqnames, strand, start and end columns are used as corresponding elements 
#' in the GRanges object. Remaining columns are added into the elementMetadata data frame.
#'
#' @param df Input data frame (required columns: seqnames, start, end, strand)
#'
#' @return GRanges object construct from df.
#' @export
dataFrameToGRanges <- function(df){
  #Convert a data.frame into a GRanges object
  
  gr = GenomicRanges::GRanges(seqnames = df$seqnames, 
                              ranges = IRanges::IRanges(start = df$start, end = df$end),
                              strand = df$strand)
  
  #Add metadata
  meta = dplyr::select(df, -start, -end, -strand, -seqnames)
  GenomicRanges::elementMetadata(gr) = meta
  
  return(gr)
}


#' Split vector of n elements into batches
#' 
#' Given a total number of elements n and batch_size, contruct a vector of length n where
# each element occurs at most batch_size times.
#'
#' @param n Length of the vector
#' @param batch_size Size of each batch
#'
#' @return vector of integers where each integer occurs at most batch_size times.
#' @export
#'
#' @examples
#' splitIntoBatches(11,3)
splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}


