import numpy as np
from tqdm.contrib import tenumerate
import pandas as pd
from functools import partial
import polars as pl
import glob
import re
import os
import argparse
from argparse import RawDescriptionHelpFormatter

# Define helper function
def get_KaKs_region(geneid_df, threshold=1, step_size=6, block_length=30):
    vec = (geneid_df['Ka/Ks'] > threshold).to_numpy()
    
    n_bases = 0
    KaKs_r = 0
    n_blocks = 0
    pos_info = "pos"
    
    if vec.sum() == 0:
        return n_bases, KaKs_r, n_blocks, pos_info
    
    n_bases = vec.sum() * step_size
    KaKs_r = vec.sum() / len(vec)
    
    for start in range(len(vec)):
        if vec[start]:
            for end in range(start, len(vec)):
                length = end - start + 1
                if vec[start:end+1].sum() != length:
                    if (end - start) >= (block_length // step_size):
                        n_blocks += 1
                        pos_info += f"_s{start * step_size}_e{(end) * step_size}"
                    vec[start:end+1] = False
                    break
    
    return n_bases, KaKs_r, n_blocks, pos_info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate kaks results from kaks calculator and find high KaKs region (>= 2 KaKs values). \n" +
        "\n"
        "Example usage: \n" +
        "python concatenate_and_find_high_KaKs.py \n" + 
        "--kaks_dir KaKs_pipeline_result_250622/strict/C1G11_C3D7v2 \n" + 
        "--sample_pair C1G11_C3D7v2 --out_dir KaKs_high_250622_csv \n" +
        "--all_kaks_csv_name split_114_6.axt.kaks.all  \n" +
        "--high_kaks_csv_name KaKs_high_C1G11_C3D7v2_250622.csv",
        formatter_class=RawDescriptionHelpFormatter
    )
    parser.add_argument("--kaks_dir", required=True, type=str, help = "Enter the path to input kaks file directory e.g. KaKs_pipeline_result_250622/strict/C1G11_C1C5v2")
    parser.add_argument("--sample_pair", required=True, type=str, help = "Enter the sample pair e.g. C1G11_C1C5v2")
    parser.add_argument("--out_dir", required=True, type=str, help = "Enter the directory to which high kaks region csv data is stored e.g. KaKs_high_250622_csv. If it does not exist, it'll create one")
    parser.add_argument("--all_kaks_csv_name", default="split_114_6.axt.kaks.all", type=str, help = "Enter the output csv file name. Default argument is split_114_6.axt.kaks.all, which results in samplepair_split114_6.axt.kaks.all.csv in final name")
    parser.add_argument("--high_kaks_csv_name", required=True, type=str, help = "Enter the output csv file name. Example argument: KaKs_high_C1G11_C1C5v2_250622.csv")
    
    args = parser.parse_args()

    kaks_dir = args.kaks_dir
    sample_pair = args.sample_pair
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    all_kaks_csv_name = args.all_kaks_csv_name
    high_kaks_csv_name = args.high_kaks_csv_name

    kaks_files = kaks_dir + "/*split_114_6.axt.kaks"
    files = glob.glob(kaks_files)

    print("start concatenating kaks files")
    
    i = 0
    df_pl_all = []
    for file in files:
        i += 1
        if i % 1000 == 0:
            print(f"concatenated {i} files")
        new_df_pl = pl.read_csv(file, separator = "\t", has_header = True, infer_schema = False)
        new_df_pl = new_df_pl.with_columns(pl.col(pl.String).replace("NA", None))
        df_pl_all.append(new_df_pl)
    df_pl = pl.concat(df_pl_all, how="vertical_relaxed")

    df_pd = df_pl.to_pandas()
    df_pd['pos'] = df_pd['Sequence'].str.extract(r'\((\d+)\-', expand=False).astype(int)
    df_pd['pair'] = df_pd['Sequence'].str.extract(r'.{11}(.*)\(', expand=False)
    df_pd['geneid'] = df_pd['Sequence'].str.extract(r'(.{10})', expand=False)
    df_pd["Ka/Ks"] = df_pd["Ka/Ks"].apply(pd.to_numeric, errors='coerce')

    all_kaks_csv_path = out_dir + "/" + sample_pair + all_kaks_csv_name + ".csv"
    df_pd.to_csv(all_kaks_csv_path, index=False)

    geneid_vec = df_pd['geneid'].unique()

    # initialize KaKs > 1 data frame
    KaKs_all_df = pd.DataFrame({'geneid': geneid_vec})
    KaKs_all_df['KaKs_bases'] = 0
    KaKs_all_df['KaKs_r'] = 0
    KaKs_all_df['n_blocks'] = 0
    KaKs_all_df['pos_info'] = ''

    # initialize KaKs > 2 data frame
    KaKs_all_df2 = pd.DataFrame({'geneid': geneid_vec})
    KaKs_all_df2['KaKs_bases2'] = 0
    KaKs_all_df2['KaKs_r2'] = 0
    KaKs_all_df2['n_blocks2'] = 0
    KaKs_all_df2['pos_info2'] = ''

    print("start calculating high KaKs region")
    for idx, geneid in tenumerate(geneid_vec):
        geneid_df = df_pd[df_pd['geneid'] == geneid]
        res1 = get_KaKs_region(geneid_df)
        KaKs_all_df.loc[KaKs_all_df['geneid'] == geneid, ['KaKs_bases', 'KaKs_r', 'n_blocks', 'pos_info']] = res1
        res2 = get_KaKs_region(geneid_df, threshold=2)
        KaKs_all_df2.loc[KaKs_all_df2['geneid'] == geneid, ['KaKs_bases2', 'KaKs_r2', 'n_blocks2', 'pos_info2']] = res2
        #print(f"finished {geneid}")

    KaKs_all_df_combined = pd.merge(KaKs_all_df, KaKs_all_df2, on='geneid')

    # Specify the output csv file name and directory
    high_kaks_csv_path = out_dir + "/" + high_kaks_csv_name
    KaKs_all_df_combined.to_csv(high_kaks_csv_path, index=False)





