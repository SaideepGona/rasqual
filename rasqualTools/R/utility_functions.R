#Utility functions that are often shared between multiple pacakages, but are too small to deserve their own package.

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

#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param A instance of GRanges, RangedData, or RangesList 
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}

