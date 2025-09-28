## Get the list of failed genes 

import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import argparse

def test_multiple_three(gap_str):
    if type(gap_str) is str:
        gap_blocks = gap_str.split(",")
        for gap_block in gap_blocks:
            #print(gap_block)
            gap_start_stop = gap_block.split("-")
            gap_len = int(gap_start_stop[1]) - int(gap_start_stop[0]) + 1
            if gap_len % 3 == 0:
                continue
            else:
                return "No"
        return "Yes"
    else:
        return "Yes"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get gap information in the fasta files in given directory"
    )
    parser.add_argument("--out_dir", required=True, type=str, 
                        help = "Enter the directory for output failed gene list")
    parser.add_argument("--out_fname", required=False, type=str, default = "failed_genes_all.txt",
                        help = "Enter the file name storing the failed cds gene information")
    parser.add_argument("--cds_failed_fpath", required=True, type=str, 
                        help = "Enter the path to the file storing failed_cds")
    parser.add_argument("--gap_info_dir", required = True, type=str,
                        help = "Enter the directory for the gap information")

    args = parser.parse_args()

    out_dir = args.out_dir
    out_fname = args.out_fname
    cds_failed_fpath = args.cds_failed_fpath
    gap_info_dir = args.gap_info_dir

    os.makedirs(out_dir, exist_ok=True)


    REF_df = pd.read_csv(os.path.join(gap_info_dir, "REF_gaps.csv"))
    C1G11_df = pd.read_csv(os.path.join(gap_info_dir, "C1G11_gaps.csv"))
    C1C5v2_df = pd.read_csv(os.path.join(gap_info_dir, "C1C5v2_gaps.csv"))
    C3D7v2_df = pd.read_csv(os.path.join(gap_info_dir, "C3D7v2_gaps.csv"))



    gap_combined = REF_df.merge(C1G11_df,
                                        left_on='gene_id', 
                                        right_on='gene_id')
    gap_combined = gap_combined.merge(C1C5v2_df,
                                        left_on='gene_id', 
                                        right_on='gene_id')
    gap_combined = gap_combined.merge(C3D7v2_df,
                                        left_on='gene_id', 
                                        right_on='gene_id')

    gap_combined["REF_multiple_of_three"] = gap_combined["REF"].apply(test_multiple_three)
    gap_combined["C1G11_multiple_of_three"] = gap_combined["C1G11"].apply(test_multiple_three)
    gap_combined["C1C5v2_multiple_of_three"] = gap_combined["C1C5v2"].apply(test_multiple_three)
    gap_combined["C3D7v2_multiple_of_three"] = gap_combined["C3D7v2"].apply(test_multiple_three)

    gap_error = gap_combined[(gap_combined["REF_multiple_of_three"] != "Yes") | (gap_combined["C1G11_multiple_of_three"] != "Yes") | (gap_combined["C1C5v2_multiple_of_three"] != "Yes") | (gap_combined["C3D7v2_multiple_of_three"] != "Yes")]
    n_error = len(gap_error["gene_id"])
    print(f"Number of genes with indels not multiple of three: {n_error}")

    with open(cds_failed_fpath, "r") as f:
        failed_cds = f.read().splitlines()
    n_error = len(failed_cds)
    print(f"Number of genes with failed cds: {n_error}")
    error_geneid_list_all = set(gap_error["gene_id"]).union(set(failed_cds))
    n_error = len(error_geneid_list_all)
    print(f"Number of genes with errors in total: {n_error}")
    with open(os.path.join(out_dir, out_fname), "w") as f:
        for gene in error_geneid_list_all:
            line = gene + "\n"
            f.write(line)