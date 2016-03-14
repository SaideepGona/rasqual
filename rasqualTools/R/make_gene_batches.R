#' Split gene ids into batches for runRasqual.py script
#'
#' @param gene_metadata Data frame with gene metadata (gene_id column required)
#' @param batch_size Number of genes in a batch
#' @param batch_prefix Prefix of the batch id. Useful if combining results from 
#' multiple calls to rasqualConstructGeneBatches
#'
#' @return Data frame of gene batches (columns: batch_id, gene_ids)
#' @export
rasqualConstructGeneBatches <- function(gene_metadata, batch_size, batch_prefix = "batch"){
  batch_df = dplyr::select(gene_metadata, gene_id) %>%
    dplyr::mutate(batch_number = splitIntoBatches(length(gene_id), batch_size)) %>%
    dplyr::mutate(batch_id = paste(batch_prefix, batch_number, sep = "_")) %>% 
    dplyr::group_by(batch_id) %>%
    dplyr::summarize(gene_ids = paste(gene_id, collapse = ","))
  return(batch_df)
}


#' Split genes into different batch sizes based on how many cis and feature SNPs they have.
#'
#' This script calculates the number of tests required for each gene (feature_snp_count * cis_snp_count)
#' and splits genes into four groups based on this value (<15,000; 15,000 to 50,000; 50,000 to 100,000 and > 100,000).
#' @param gene_metadata Data frame with gene meta data, required columns: gene_id, feature_snp_count, cis_snp_count.
#' @param batch_sizes Vector of length 4. Contains the number of gene in a batch for each of the four groups.
#' @param batch_prefix Prefix of the batch id, default "batch".
#'
#' @return Data frame of genes split into batches.
#' @export
rasqualOptimisedGeneBatches <- function(gene_metadata, batch_sizes = c(20,8,3,1), batch_prefix = "batch"){
  #Calculate the number of tests
  gene_metadata = dplyr::mutate(gene_metadata, test_count = feature_snp_count*cis_snp_count)
  
  #Split genes into batches based on the number of tests
  quick_genes = dplyr::filter(gene_metadata, test_count < 15000) %>% 
    rasqualConstructGeneBatches(batch_sizes[1], paste(batch_prefix, "1", sep = "_"))
  medium_genes = dplyr::filter(gene_metadata, test_count >= 15000, test_count < 50000) %>%
    rasqualConstructGeneBatches(batch_sizes[2], paste(batch_prefix, "2", sep = "_"))
  slow_genes = dplyr::filter(gene_metadata, test_count >= 50000, test_count < 100000)%>%
    rasqualConstructGeneBatches(batch_sizes[3], paste(batch_prefix, "3", sep = "_"))
  extra_slow_genes = dplyr::filter(gene_metadata, test_count >= 100000)%>%
    rasqualConstructGeneBatches(batch_sizes[4],paste(batch_prefix, "4", sep = "_"))
  
  batches = rbind(quick_genes, medium_genes, slow_genes, extra_slow_genes)
  return(batches)
}