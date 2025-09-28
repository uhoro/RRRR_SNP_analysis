

import os
import sys
import pandas as pd
from itertools import combinations
import glob
import argparse

def count_ele(vec1, vec2):
    return sum(1 for x in vec1 if x in vec2)

def get_overlap(s1, s2, data):
    is_s1_s2_vec = []
    u_s1_s2_vec = []
    u_s1_s2_f_vec = []
    u_s2_s1_vec = []
    u_s2_s1_f_vec = []
    
    for i, row in data.iterrows():
        #if i % 100000 == 0:
            #print(f"finished processing variant: {i}")
        vec1 = row[s1].split('/')
        #print("vec1", vec1)
        vec2 = row[s2].split('/')
        #print("vec2", vec2)
        
        is_s1_s2 = list(set(vec1) & set(vec2))
        is_s1_s2_vec.append(''.join(is_s1_s2))
        #print("is_s1_s2_vec", is_s1_s2_vec)
        
        u_s1_s2 = list(set(vec1) - set(vec2)) # get allele that is unique to s1
        u_s1_s2_vec.append(''.join(u_s1_s2))
        #print("u_s1_s2_vec", u_s1_s2_vec)
        
        u_s1_s2_f = count_ele(vec1, u_s1_s2)
        u_s1_s2_f_vec.append(u_s1_s2_f)
        #print("u_s1_s2_f_vec", u_s1_s2_f_vec)
        
        u_s2_s1 = list(set(vec2) - set(vec1))
        u_s2_s1_vec.append(''.join(u_s2_s1))
        #print("u_s2_s1_vec", u_s2_s1_vec)
        
        u_s2_s1_f = count_ele(vec2, u_s2_s1)
        u_s2_s1_f_vec.append(u_s2_s1_f)
        #print("u_s2_s1_f_vec", u_s2_s1_f_vec)
    
    df = pd.DataFrame({
        f"is_{s1}_{s2}": is_s1_s2_vec,
        f"u_{s1}_{s2}": u_s1_s2_vec,
        f"u_{s1}_{s2}_f": u_s1_s2_f_vec,
        f"u_{s2}_{s1}": u_s2_s1_vec,
        f"u_{s2}_{s1}_f": u_s2_s1_f_vec
    })
    return df


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get the variant overlap for a given pair of samples"
    )
    parser.add_argument("--data_fname", required=True, type=str, help = "Enter path to the variant data file ")
    parser.add_argument("--s1", required=True, type=str, help = "Enter the sample1 id. " \
                                        "E.g. SM44B")
    parser.add_argument("--s2", required=True, type=str, help = "Enter the sample2 id." \
                                        "E.g. C1C5v2")
    parser.add_argument("--pair_overlap_dir", required=True, type=str, 
                        help = "Enter the dir where result dataframe is stored." \
                                "bindedp4_set240924_filt250821_INDEL_pair_overlap")
    args = parser.parse_args()
    data_fname = args.data_fname
    s1= args.s1
    s2 = args.s2
    pair_overlap_dir = args.pair_overlap_dir

    os.makedirs(pair_overlap_dir, exist_ok=True)

    # log_fname = s1 + "_" + s2 + "_pair_log"
    # file = open(os.path.join(pair_overlap_dir, log_fname), "w")
    # sys.stderr.write = file.write
    
    data = pd.read_csv(data_fname)
    for colname in list(data):
        data[colname] = data[colname].astype(str)

    result = get_overlap(s1, s2, data)
    
    result_fname = s1 + "_" + s2 + "_overlap.csv"
    result.to_csv(os.path.join(pair_overlap_dir, result_fname), 
                  index=False)

