
### Installation
The easiest way to install rasqualTools is to use devtools:

```r
install.packages("devtools")
library("devtools")
devtools::install_github("kauralasoo/rasqual/rasqualTools")
```

rasqualTools imports the following other packages: Rsamtools, readr, dplyr, GenomicRanges, plyr, magrittr, assertthat.

### STEP 1: Save read count matrix onto disk
Our read count table is on the `counts_matrix` object:

```r
print(counts_matrix[1:6, 1:6])
```

```
##                 sample_1 sample_2 sample_3 sample_4 sample_5 sample_6
## ENSG00000152213     1810      942      119      116      763      644
## ENSG00000154642      165      435      777     1895      236      874
## ENSG00000101825      348      161      290     1602      125      215
## ENSG00000089057    10289     7205     1679     1829     8934     8137
## ENSG00000272842        3        0        0        0        2        0
## ENSG00000254201        0        0        0        2        0        0
```
RASQUAL needs the counts matrix to be in a binary format. We can use `saveRasqualMatrices` command to save poth plain text as well as binary versions of the matrixt to disk:

```r
saveRasqualMatrices(list(cellTypeA = counts_matrix), "../output/", file_suffix = "expression")
```

```
## [1] "../output//cellTypeA.expression.txt"
```

### STEP 2: Calculate size factors
The simples option is to just use the library size of each sample as the size factor:

```r
size_factors = rasqualCalculateSampleOffsets(counts_matrix, gc_correct = FALSE)
saveRasqualMatrices(list(cellTypeA = size_factors), "../output/", file_suffix = "size_factors")
```

```
## [1] "../output//cellTypeA.size_factors.txt"
```
Alternatively, if we have a `gene_metadata` data frame that matches gene_ids to their average gc content then we can use that information to correct for the differences in GC bias between samples.

```r
dplyr::select(gene_metadata, gene_id, percentage_gc_content)
```

```
## Source: local data frame [2,000 x 2]
## 
##            gene_id percentage_gc_content
##              (chr)                 (dbl)
## 1  ENSG00000152213                 47.31
## 2  ENSG00000154642                 37.16
## 3  ENSG00000101825                 43.18
## 4  ENSG00000089057                 44.93
## 5  ENSG00000272842                 26.79
## 6  ENSG00000254201                 41.72
## 7  ENSG00000253557                 40.76
## 8  ENSG00000100014                 42.29
## 9  ENSG00000100836                 50.48
## 10 ENSG00000169100                 58.91
## ..             ...                   ...
```

```r
size_factors = rasqualCalculateSampleOffsets(counts_matrix, gene_metadata, gc_correct = TRUE)
saveRasqualMatrices(list(cellTypeA = size_factors), "../output/", file_suffix = "size_factors_gc")
```

```
## [1] "../output//cellTypeA.size_factors_gc.txt"
```

### STEP 3: Calculate the number of SNPs overlapping each gene
We can use the `countSnpsOverlapingExons` function to do that. Lets look at a small example. Here is the gene_metaata data frame with minimal required columns:

```r
gene_data = dplyr::select(gene_metadata, gene_id, chr, strand, exon_starts, exon_ends)[c(1,6,20,34),] 
print(gene_data)
```

```
## Source: local data frame [4 x 5]
## 
##           gene_id   chr strand       exon_starts         exon_ends
##             (chr) (chr)  (int)             (chr)             (chr)
## 1 ENSG00000152213    13      1 49628299,49630430 49628615,49633872
## 2 ENSG00000254201     8      1 19246350,19249163 19246450,19249240
## 3 ENSG00000183186    19     -1     405438,409006     408401,409139
## 4 ENSG00000198816    19      1   7516118,7519202   7516249,7521026
```
And here are the SNP coordinates:

```r
print(snp_coords)
```

```
## Source: local data frame [8,362 x 3]
## 
##      chr    pos      snp_id
##    (chr)  (dbl)       (chr)
## 1     19 226776  rs76534612
## 2     19 230130 rs200141179
## 3     19 240867   rs1975526
## 4     19 244421  rs62103026
## 5     19 244426  rs62103027
## 6     19 245410 rs140409643
## 7     19 245631  rs67286684
## 8     19 245680 rs148381057
## 9     19 245844  rs59441037
## 10    19 245867  rs57114567
## ..   ...    ...         ...
```
We can now count how many SNP fall within the gene itself or its cis window. This command adds the `feature_snp_count` and `cis_snp_count` columns to the imput gene_metadata data frame.

```r
snp_counts = countSnpsOverlapingExons(gene_data, snp_coords, cis_window = 5e5)
dplyr::select(snp_counts, gene_id, feature_snp_count, cis_snp_count)
```

```
## Source: local data frame [4 x 3]
## 
##           gene_id feature_snp_count cis_snp_count
##             (chr)             (int)         (int)
## 1 ENSG00000152213                 0             0
## 2 ENSG00000183186                 8          1881
## 3 ENSG00000198816                 2          2850
## 4 ENSG00000254201                 0             0
```

