

from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import pandas as pd
import os
import glob
import argparse
import shutil


def fasta2axt_pairs(fasta_file, sample_ids, out_dir):

    file_name = str.split(fasta_file, sep = "/")[-1]
    sample_gene_id = str.split(file_name, sep = ".")[0]
    gene_num = str.split(sample_gene_id, sep = "_")[-1]
    gene_id = "FUN_" + gene_num

    sequence = {}
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            sequence[record.id] = str(record.seq)

    for i in range(len(sample_ids)):
        for j in range(i+1, len(sample_ids)):
            sample1_id = sample_ids[i]
            sample2_id = sample_ids[j]
            os.makedirs(os.path.join(out_dir, sample1_id + "_" + sample2_id), exist_ok=True)

            sample1_seq = sequence[sample1_id]
            sample2_seq = sequence[sample2_id]
            sinkfname = os.path.join(out_dir, sample1_id  + "_" + sample2_id, \
                                        sample1_id + "_" + sample2_id + "_" + gene_id + ".axt")
            with open(sinkfname, 'w') as f:
                pass  # Clear the file

            # Write results to file
            with open(sinkfname, 'a') as f:
                header = f"{gene_id}_{sample1_id}_{sample2_id}"
                f.write(header + "\n")
                f.write(str(sample1_seq) + "\n")
                f.write(str(sample2_seq) + "\n")
                f.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get axt file from fasta file"
    )
    parser.add_argument("--initial_fasta_dir", required=True, type=str, help = "Enter the directory where the fasta sequences are stored. " \
                                        "E.g. out/fasta_files_strict_250823")
    parser.add_argument("--clustalo_fasta_dir", required=True, type=str, help = "Enter the directory where the clustalo alignmed fasta sequences are stored." \
                                        "E. g. clustalo/result/fasta_files_strict_250823")
    # parser.add_argument("--genes_with_errors_file", required=True, type=str, help = "Enter the text file name where the genes with errors after clustalo is stored." \
    #                                     "E.g. out/mutations_df/strict_250823_after_clustalo_genes_with_errors.txt")
    parser.add_argument("--curated_fasta_dir", required=True, type=str, help = "Enter output dir name for curated fasta files" \
                        "E.g. out/fasta_files_strict_250823_curated")
    parser.add_argument("--out_axt_dir", required=True, type=str, help = "Enter output dir name for axt files" \
                        "E.g. out/axt_files_strict_250823")
    parser.add_argument("--sample_list", nargs='+', default = ["C1G11", "C1C5v2", "C3D7v2"],
                        help = "Enter sample to be considered from which star is removed")
    parser.add_argument("--test", required=False, type=bool, default = False, help = "Enable test mode with True")
    args = parser.parse_args()


    initial_fasta_file_dir = args.initial_fasta_dir
    initial_fasta_file_dir += "/*"
    initial_files = glob.glob(initial_fasta_file_dir)
    clustalo_fasta_dir = args.clustalo_fasta_dir
    clustalo_fasta_dir_w = clustalo_fasta_dir + "/*"
    clustalo_files = glob.glob(clustalo_fasta_dir_w)
    #genes_with_errors_file = args.genes_with_errors_file
    curated_fasta_dir = args.curated_fasta_dir
    out_axt_dir = args.out_axt_dir
    sample_list = args.sample_list
    test = args.test

    os.makedirs(curated_fasta_dir, exist_ok=True)

    #with open(genes_with_errors_file, 'r') as file:
    #    remaining_error_genes = [line.strip() for line in file]

    clustalo_filenames = os.listdir(clustalo_fasta_dir)
    clustalo_heads = [f.split(".")[0] for f in clustalo_filenames]
    clustalo_genes = ["FUN_" + file_head.split("_")[-1] for file_head in clustalo_heads]
    print(f"Clustalo_genes_example: {clustalo_genes[0]}")

    # Prepare curated_fasta_dir replacing error gene fasta file with clustalo aligned fasta file
    i = 0
    initial_files.sort()
    for file in initial_files:
        i += 1
        if test:
            if i > 50:
                break
        if i % 1000 == 0:
            print(f"processed {i} th gene")
        file_name = str.split(file, "/")[-1]
        print(f"file_name: {file_name}")
        file_head = str.split(file_name, ".")[0]
        print(f"file_head: {file_head}")
        gene_id = "FUN_" + file_head.split("_")[-1]
        if gene_id in clustalo_genes:
            print(f"include clustalo {gene_id}")
            clustalo_gene_file = [f for f in clustalo_files if gene_id in f][0]
            print(clustalo_gene_file)
            new_file_name = str.split(clustalo_gene_file, sep = "/")[-1]
            print(new_file_name)
            new_file_name = str.split(new_file_name, sep = ".")[0] + ".fa"
            new_file_path = os.path.join(curated_fasta_dir, new_file_name)
            shutil.copy(clustalo_gene_file, new_file_path)
        else:
            new_file_path = os.path.join(curated_fasta_dir, file_name)
            shutil.copy(file, new_file_path)

    curated_fasta_dir_w = curated_fasta_dir + "/*"
    fasta_files = glob.glob(curated_fasta_dir_w)
    i = 0
    for fasta_file in fasta_files:
        i += 1
        if test:
            if i > 10:
                break
        if i % 1000 == 0:
            print(f"start generating axt files for ${i} th gene")
        fasta2axt_pairs(fasta_file, sample_list, out_axt_dir)


