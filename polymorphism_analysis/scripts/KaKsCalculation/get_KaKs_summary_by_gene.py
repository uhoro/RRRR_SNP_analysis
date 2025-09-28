import numpy as np
from tqdm.contrib import tenumerate
import pandas as pd
from functools import partial
import polars as pl
import glob
import re
import os
import argparse
from os import listdir
from os.path import isfile, join

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="get KaKs summary by gene from KaKs all concatenated csv files (output of concatenate_kaks_results.py)"
    )
    #parser.add_argument("--gene_list", required=True, action="append", help = "Enter the list of genes whose summary to be obtained, e.g. FUN_001390, FUN_001442, FUN_002102")
    parser.add_argument("--gene_list", required=True, nargs='+', help = "Enter the list of genes whose summary to be obtained, e.g. FUN_001390 FUN_001442 FUN_002102")
    parser.add_argument("--sample_pairs", nargs='+', default = ["C1G11_C1C5v2", "C1G11_C3D7v2", "C1C5v2_C3D7v2"],  
                        help = "Enter the sample pair list. By default, C1G11_C1C5v2 C1G11_C3D7v2 C1C5v2_C3D7v2")
    parser.add_argument("--kaks_all_dir", required=True, type=str, help = "Enter the directory where all kaks data for sample pair is stored.")
    parser.add_argument("--kaks_all_file_suffix", default = "split_114_6.axt.kaks.all", type=str, help = "Enter the kaks all csv file name suffix. By default, split_114_6.axt.kaks.all")
    parser.add_argument("--out_folder_path", default = "KaKs_pipeline_result/kaks_summary_by_gene", type=str, help = "Enter the path to output folder. By default, KaKs_pipeline_result/kaks_summary_by_gene ")
    args = parser.parse_args()

    gene_list = args.gene_list
    sample_pairs = args.sample_pairs
    kaks_all_dir = args.kaks_all_dir
    kaks_all_file_suffix = args.kaks_all_file_suffix
    out_folder_path = args.out_folder_path
    os.makedirs(out_folder_path, exist_ok=True)
    print(f"out folder path: {out_folder_path}")

    df_dict = dict.fromkeys(sample_pairs)
    for sample_pair in sample_pairs:
            print("start processing ", sample_pair)
            csv_name = kaks_all_dir + "/" + sample_pair + kaks_all_file_suffix + ".csv"
            df_dict[sample_pair] = pd.read_csv(csv_name)

    for gene_of_interest in gene_list:
        print("start processing ", gene_of_interest)
        df_list = []
        for sample_pair in sample_pairs:
            print("start processing ", sample_pair)
            samplepair_df_all = df_dict[sample_pair]
            samplepair_df_gene = samplepair_df_all[samplepair_df_all["geneid"] == gene_of_interest]
            print("size of df: ", str(samplepair_df_gene.shape))
            df_list.append(samplepair_df_gene)
        df_list_all = pd.concat(df_list, ignore_index = True)
        df_list_agg = df_list_all.pivot(index = "pos", columns = "pair", values = "Ka/Ks")
        df_list_agg.columns = [col + "_KaKs" for col in df_list_agg.columns]
    
        df_list_agg["avg_KaKs"] = df_list_agg.mean(axis = 1)
        df_list_agg = df_list_agg.reset_index()
        print("size of final df: ", str(df_list_agg.shape))
        out_fname = out_folder_path + "/" + gene_of_interest + ".csv"
        print("outfname: ", out_fname)
        df_list_agg.to_csv(out_fname, index = False)

