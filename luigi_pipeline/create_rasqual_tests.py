'''
Author : Saideep Gona

Imports a windowed, intersected GTF,VCF file and pulls out the relevant components for 
performing rasqual tests.
'''

import sys, os, argparse
import pickle

cwd = os.getcwd()
print("Python Interpreter: ", sys.executable)

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Path to windowed, intersected, GTF-VCF file")
parser.add_argument("output_file", help="Path to file of rasqual tests")
parser.add_argument("raw_count_table", help="Path to file raw counts. The features must match the rows in this table")
args = parser.parse_args()

# INTERSECT LINE
# '''
# 1
# transcribed_unprocessed_pseudogene
# transcript
# 10872
# 15412
# .
# +
# .
# gene_id "ENSG00000223972"; transcript_id "ENST00000515242"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-201"; transcript_source "ensembl";
# 1  (9)
# 13380
# 1:13380:C:G	C
# G	
# .
# .
# AC=2;AN=22
# GT
# 1|1	0|0	0|0	0|00|0	0|0	0|0	0|0	0|0	0|0	0|0
# '''

if not args.raw_count_table.endswith(".csv"):
    print("Raw count table should be a .csv file")
    sys.exit()

valid_feature_set = set()
valid_features = []
count = 0
with open(args.raw_count_table, "r") as raw:
    for line in raw:
        if count == 0:
            count += 1
            continue
        p_line = line.rstrip("\n").split(",")
        valid_feature_set.add(p_line[0])
        valid_features.append(p_line[0])



all_feature_dict = {}
all_features = []
n_samples = 0

with open(args.input_file, "r") as inp:
    gene_count = 0
    for line in inp:
        # print(line)
        p_line = line.rstrip("\n").split("\t")

        n_samples = str(len(p_line[-1].split(" ")))
        
        loci = " ".join([p_line[0],p_line[3],p_line[4]]) 

        snp = "_".join([p_line[9] + p_line[10]])

        gtf_params = p_line[8].split("; ")
        gtf_params_dict = {} 
        for param in gtf_params:
            spl = param.split(" ")
            gtf_params_dict[spl[0]] = spl[1].strip('"')
        gene_name = gtf_params_dict["gene_name"]


        if p_line[2] == "gene":
            if gene_name not in all_feature_dict:
                gene_count += 1
                print("Gene Count: ", gene_count, " ", gtf_params_dict["gene_name"])
                all_features.append(gene_name)
                all_feature_dict[gene_name] = {
                    "chrom": p_line[0],
                    "start": p_line[3],
                    "end": p_line[4],
                    "snps": {snp},
                    "feature_snps": set(),
                    "exons": set(),
                    "exon_starts": [],
                    "exon_ends": [],
                    "gene_id": gtf_params_dict["gene_id"]
                }
            else:
                all_feature_dict[gene_name]["snps"].add(snp)

        if p_line[2] == "exon":
            all_feature_dict[gene_name]["feature_snps"].add(snp)
            if loci not in all_feature_dict[gene_name]["exons"]:
                all_feature_dict[gene_name]["exons"].add(loci)
                all_feature_dict[gene_name]["exon_starts"].append(loci.split(" ")[1])
                all_feature_dict[gene_name]["exon_ends"].append(loci.split(" ")[2])

feature_count = 1

try:

    with open(args.output_file, "w") as out:

        for feature in valid_features:
            if feature not in all_feature_dict:
                raise RuntimeError("Features from GTF and count matrix do not match")
                
            loc = [
                all_feature_dict[feature]["chrom"],
                ":",
                all_feature_dict[feature]["start"],
                "-",
                all_feature_dict[feature]["end"]
            ]
            rasqual_command = [
                "PLACEHOLDER_TABIX", "PLACEHOLDER_VCF",
                "".join(loc),
                "|",
                "PLACEHOLDER_RASQUAL",
                "-y", "PLACEHOLDER_YBIN",
                "-k", "PLACEHOLDER_KBIN",
                "-n", n_samples,
                "-j", str(feature_count),
                "-l", str(len(all_feature_dict[feature]["snps"])),
                "-m", str(len(all_feature_dict[feature]["feature_snps"])),
                "-s", ",".join(all_feature_dict[feature]["exon_starts"]),
                "-e", ",".join(all_feature_dict[feature]["exon_ends"]),
                "-f", feature,
                "--n-threads", "PLACEHOLDER_NUMTHREADS"
            ]

            feature_count += 1

            out.write(" ".join(rasqual_command) + "\n")

except:
    os.remove(args.output_file)

print("FINISHED COLLECTING SNP/FEATURE INFO")




        
        