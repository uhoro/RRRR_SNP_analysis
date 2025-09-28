import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
import os
import argparse
from os import listdir
from os.path import isfile, join

def get_gap_info(filepath, sampleid):
    """_summary_

    Args:
        filepath (_type_): _description_
            e.g. fastapath = "/Users/uhoro/bio/pasteur/variant_calling/notebooks/250828_analysis/data/strict_250823_curated_FUN_007197.fa"
        sampleid (_type_): _description_
            e.g. "REF"
    """
    
    ref = get_sample_sequence_fast(filepath, sampleid)
    seq_str = str(ref.seq)
    #seq_str = "ATG--CAT---ACGT"
    seq_str_lis = list(seq_str)

    in_gap = False
    gap_start = 0
    gap_end = 0
    gap_res = []
    for i in range(len(seq_str_lis)):
        if in_gap:
            if seq_str_lis[i] == "-":
                continue
            else:
                gap_end = i
                gap_info = str(gap_start) + "-" + str(gap_end)
                gap_res.append(gap_info)
                in_gap = False
        else:
            if seq_str_lis[i] == "-":
                in_gap = True
                gap_start = i + 1
    gap_res_str = ",".join(gap_res)
    return gap_res_str


def get_sample_sequence_fast(filepath, sampleid):
    with open(filepath, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if sampleid in record.id:
                return SeqIO.SeqRecord(Seq(str(record.seq)), id=sampleid)
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get gap information in the fasta files in given directory"
    )
    parser.add_argument("--out_dir", required=True, type=str, 
                        help = "Enter the directory where the gap information data is stored")
    parser.add_argument("--sampleid", required=True, type=str, 
                        help = "Enter sample id for which gap information is retrieved")
    # parser.add_argument("--sample_list", nargs='+', default = ["C1G11", "C1C5v2", "C3D7v2"],
    #                     help = "Enter sample list to be selected. Default is C1G11 C1C5v2 C3D7v2")
    parser.add_argument("--data_dir", required = True, type=str,
                        help = "Enter the directory of original fasta file")
    parser.add_argument("--test", required=False, type=bool, default = False, 
                        help = "Enable test mode with True")
    args = parser.parse_args()

    out_dir = args.out_dir
    sampleid = args.sampleid
    data_dir = args.data_dir
    test = args.test

    os.makedirs(out_dir, exist_ok=True)

    original_files = [f for f in listdir(data_dir) if isfile(join(data_dir, f))]
    print(f"Number of seqs in {data_dir}: {len(original_files)}")
    # for test, do only for the first 10 geneid_vec
    if test:
        last = 10
    else:
        last = len(original_files)

    i = 0
    gap_info_list = [pd.DataFrame(columns=['gene_id', sampleid])]
    n_gap_genes = 0
    for file in original_files:
        if i >= last:
            break
        i += 1
        if i % 1000 == 0:
            print(f"{i} th gene")
            print(f"start processing: {file}")
        geneid = "FUN_" + file.split(".")[0].split("_")[-1]
        gap_info = get_gap_info(data_dir + "/" + file, sampleid)
        new_row = pd.DataFrame({'gene_id': [geneid], 
                                sampleid : [gap_info]})
        gap_info_list.append(new_row)
        if gap_info != "":
               n_gap_genes += 1
  
    df_gap = pd.concat(gap_info_list, ignore_index = True)
    out_fname = sampleid + "_gaps.csv"
    df_gap.to_csv(os.path.join(out_dir, out_fname), index = False)
    print(f"Number of genes considered: {len(original_files)}")
    print(f"Number of gap genes: {n_gap_genes}")


