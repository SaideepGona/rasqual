'''
Author : Saideep Gona

This script is designed to modify the final gene counts 
matrix such that the rownames are syntactically correct
for downstream analysis.
'''

import sys, os, argparse

cwd = os.getcwd()
print("Python Interpreter: ", sys.executable)

parser = argparse.ArgumentParser()
parser.add_argument("input_path", help="Input path")
parser.add_argument("output_path", help="Output path")
args = parser.parse_args()

with open(gene_count_ori, "r") as gco:
    with open(assay_out, "w") as ao:

        line_num = 0
        for line in gco:
            if line_num == 0:
                line_num += 1
                header = line.rstrip("\n").split(",")
                new_header = [x.split("_")[0] for x in header]
                ao.write(",".join(new_header) + "\n")
                continue
            
            p_line = line.rstrip("\n").split(",")
            genes = p_line[0]
            if len(genes.split("|")) == 2:
                p_line[0] = genes.split("|")[-1]
            elif len(genes.split("|")) > 2:
                p_line[0] = genes.split("|")[0]

            ao.write(",".join(p_line) + "\n")

