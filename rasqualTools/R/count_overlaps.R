#' Count the number of feature SNPs and cis SNPs overlapping a set of peak calls
#'
#' @param peak_metadata Data frame with peak metadata (required columns: gene_id, chr, start, end)
#' @param snp_coords Data frame with SNP coordinates from a VCF file (required columns: chr, pos, snp_id)
#' @param cis_window Size of the cis window from both sides of the peak.
#'
#' @return Data frame with peak coordinates, cis region coordiantes as well as number of cis and feature snps.
#' @export
countSnpsOverlapingPeaks <- function(peak_metadata, snp_coords, cis_window = 500000){
  
  #Construct peak coords data frame
  peak_coords = dplyr::transmute(peak_metadata, gene_id, chromosome_name = chr, strand = 1, 
                                 exon_starts = start, exon_ends = end) %>%
    dplyr::mutate(range_start = pmax(0, exon_starts - cis_window), range_end = exon_ends + cis_window)
  
  #Construct GRanges objects
  peak_granges = GenomicRanges::GRanges(seqnames = peak_coords$chromosome_name, 
                                        ranges = IRanges::IRanges(start = peak_coords$exon_starts, end = peak_coords$exon_ends))
  region_granges = GenomicRanges::GRanges(seqnames = peak_coords$chromosome_name, 
                                          ranges = IRanges::IRanges(start = peak_coords$range_start, end = peak_coords$range_end))
  snp_granges = GenomicRanges::GRanges(seqnames = snp_coords$chr, ranges = IRanges::IRanges(start = snp_coords$pos, end = snp_coords$pos))
  
  #Count overlaps
  feature_snps = GenomicRanges::countOverlaps(peak_granges, snp_granges, ignore.strand = TRUE)
  cis_snps = GenomicRanges::countOverlaps(region_granges, snp_granges, ignore.strand = TRUE)
  
  new_peak_coords = dplyr::mutate(peak_coords, feature_snp_count = feature_snps, cis_snp_count = cis_snps)
  return(new_peak_coords)
}

#' Count the number of feature SNPs and cis SNPs overlapping exons of all genes
#'
#' @param gene_metadata Data frame with gene metadata (required columns: gene_id, chr, strand, exon_starts, exon_ends)
#' -exon_stars: comma-separated list of exon start coordinates
#' -exon_ends: comma-separated list of exon end coordinates.
#' @param snp_coords Data frame with SNP coordinates from a VCF file (required columns: chr, pos, snp_id)
#' @param cis_window Size of the cis window from both sides of the gene.
#'
#' @return Data frame with exon coordinates, cis region coordiantes as well as the number of cis and feature snps.
#' @export
countSnpsOverlapingExons <- function(gene_metadata, snp_coords, cis_window = 5e5){
  
  #Split exon coordinates into separate rows
  gene_df_list = plyr::dlply(gene_metadata, .(gene_id), function(x){
    data.frame(gene_id = x$gene_id, 
               seqnames = x$chr,
               strand = x$strand,
               start = as.numeric(unlist(strsplit(x$exon_starts,","))),
               end = as.numeric(unlist(strsplit(x$exon_ends,","))) )
  })
  exon_df = plyr::ldply(gene_df_list, .id = NULL) %>% tbl_df()
  
  #Counts the number of feature SNPs
  exon_granges = dataFrameToGRanges(exon_df)
  snp_granges = GenomicRanges::GRanges(seqnames = snp_coords$chr, 
                                       ranges = IRanges::IRanges(start = snp_coords$pos, end = snp_coords$pos))
  n_feature_snps = GenomicRanges::countOverlaps(exon_granges, snp_granges, ignore.strand=TRUE)
  feature_snp_df = dplyr::mutate(exon_df, feature_snp_count = n_feature_snps) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarise(seqnames = seqnames[1], strand = strand[1], start = min(start), end = max(end), 
                     feature_snp_count = sum(feature_snp_count))
  
  #Count the number of cis SNPs
  cis_df = dplyr::mutate(feature_snp_df, start = pmax(0, start - cis_window), end = end + cis_window)
  cis_granges = dataFrameToGRanges(cis_df)
  n_cis_snps = GenomicRanges::countOverlaps(cis_granges, snp_granges, ignore.strand=TRUE)
  result = dplyr::mutate(cis_df, cis_snp_count = n_cis_snps, gene_id = as.character(gene_id)) %>%
    dplyr::rename(range_start = start, range_end = end)
  
  #Add exon start and end coords and reorder columns
  start_end_df = dplyr::select(gene_metadata, gene_id, exon_starts, exon_ends)
  result = dplyr::left_join(result, start_end_df, by = "gene_id") %>%
    dplyr::transmute(gene_id, chromosome_name = seqnames, strand, exon_starts, exon_ends, 
                     range_start, range_end, feature_snp_count, cis_snp_count)
  return(result)
}