## Helper scripts
This folder contains some helper scripts that should simplify running RASQUAL on multiple samples:
* `bamCountASE.py` - Uses ASEReadCounter from the Genome Analysis Toolkit (GATK) to quantify allele-specifc expression in RNA-Seq BAM files.
* `mergeASECounts.py` - Merges allele-specifc read counts from the previous script into a single table.
* `vcfAddASE.py` - Adds allele-specifc counts from `mergeASECounts.py` into the VCF file.
* `runRasqual.py` - Run RASQUAL on a batch of genes.
* `mergeRasqualBatches.py` - Merges RASQUAL output files from multiple batches into a single text file.
* `rasqualToEigenMT.py` - Converts RASQUAL output into format suitable for [eigenMT](https://github.com/joed3/eigenMT).


