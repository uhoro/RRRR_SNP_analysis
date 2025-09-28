from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import pandas as pd
import os
import glob
import argparse
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get the dataframe of mutations for each gene. Also count the number of genes with imperfect alignment"
    )
    parser.add_argument("--fasta_dir", required=True, type=str, help = "Enter the directory where the fasta sequences are stored. " \
                                        "E.g. WD/fasta_files_JG_SW_250902_curated")
    parser.add_argument("--geneid_list_txt", required=True, type=str, help = "Enter the txt file containing list of gene id for which the number of non-synonymous mutations would be counted." \
                                        "E.g. ")
    parser.add_argument("--out_dir", required=True, type=str, help = "Enter the output directory where count matrix is saved")
    parser.add_argument("--test", required=False, type=bool, default = False, help = "Enable test mode with True")
    
    args = parser.parse_args()

    fasta_dir = args.fasta_dir
    geneid_list_txt = args.geneid_list_txt
    out_dir = args.out_dir
    test = args.test

    os.makedirs(out_dir, exist_ok=True)


    with open(geneid_list_txt) as f:
        geneid_list = f.read().splitlines() 
    n_processed = 0
    for gene_id in geneid_list:
        n_processed += 1
        if (n_processed > 10) and test:
                break
        if n_processed % 50 == 0:
            print(f"Start processing: {n_processed} th gene")
        
        file_list = os.listdir(fasta_dir)
        filename = [file for file in file_list if gene_id in file][0]
        fasta_file_full = os.path.join(fasta_dir, filename)
        # Read sequences from FASTA file
        sequences = {}
        # print(f"start processing {fasta_file}")
        for record in SeqIO.parse(fasta_file_full, "fasta"):
            sequences[record.id] = str(record.seq)

        result_table = pd.DataFrame(columns=sequences.keys(), index = sequences.keys())
        sample_list = list(sequences.keys())
        for sample1 in sample_list:
            for sample2 in sample_list:
                non_synonymous_sub = 0
                for base_i in range(0, len(sequences[sample1]), 3):
                    if ("-" in sequences[sample1][base_i:base_i+3]) or ("-" in sequences[sample2][base_i:base_i+3]):
                        continue
                    if len(sequences[sample1][base_i:base_i+3]) % 3 != 0:
                        print(f"geneid {gene_id}, base_i: {base_i}, codon: {sequences[sample1][base_i:base_i+3]}")
                    aa1 = Seq(sequences[sample1][base_i:base_i+3]).translate()
                    aa2 = Seq(sequences[sample2][base_i:base_i+3]).translate()
                    if aa1 != aa2:
                        non_synonymous_sub += 1
                result_table.loc[sample1, sample2] = non_synonymous_sub
        
        out_csv_fname = gene_id + "_NonSyn_count.csv"
        out_csv_file = os.path.join(out_dir, out_csv_fname)
        result_table.to_csv(out_csv_file, index = False)
    

