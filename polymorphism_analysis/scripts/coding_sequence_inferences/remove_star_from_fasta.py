import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
import os
import argparse
from os import listdir
from os.path import isfile, join

def remove_star(filepath, sample_lis):
    """_summary_

    Args:
        filepath (_type_): _description_
            e.g. fastapath = "/Users/uhoro/bio/pasteur/variant_calling/notebooks/250828_analysis/data/strict_250823_curated_FUN_007197.fa"
        sampleid (_type_): _description_
            e.g. "REF"
    """

    seq_trimmed_dict = {}
    for sampleid in sample_lis:
        ref = get_sample_sequence_fast(filepath, sampleid)
        seq_str = str(ref.seq)  
        seq_trimmed_dict[sampleid] = seq_str.replace("*", "")
    
    return seq_trimmed_dict


# sample_lis = ["S1", "S2"]
# seq_dict = {}
# seq_trimmed_dict = {}
# seq_dict["S1"] = "AT---GC-DFED"
# seq_dict["S2"] = "AT--AGC--GCT"
# seq_trimmed_dict["S1"] = ""
# seq_trimmed_dict["S2"] = ""
# seq_trimmed_dict

def get_sample_sequence_fast(filepath, sampleid):
    with open(filepath, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if sampleid in record.id:
                return SeqIO.SeqRecord(Seq(str(record.seq)), id=sampleid)
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove star from fasta files in given directory"
    )
    parser.add_argument("--out_dir", required=True, type=str, 
                        help = "Enter the directory where the gap removed fasta file is stored")
    parser.add_argument("--sample_list", nargs='+', default = ["REF", "C1G11", "C1C5v2"],
                        help = "Enter sample to be considered from which star is removed")
    parser.add_argument("--data_dir", required = True, type=str,
                        help = "Enter the directory of original fasta file")
    parser.add_argument("--test", required=False, type=bool, default = False, 
                        help = "Enable test mode with True")
    args = parser.parse_args()

    out_dir = args.out_dir
    sample_list = args.sample_list
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

    for file in original_files:
        if i >= last:
            break
        i += 1
        if i % 500 == 0:
            print(f"{i} th gene")
            print(f"start processing: {file}")
        geneid = "FUN_" + file.split(".")[0].split("_")[-1]
        trimmed_seq_dict = remove_star(data_dir + "/" + file, sample_list)

        out_seq = []
        for key in trimmed_seq_dict.keys():
            seq_str = trimmed_seq_dict[key]
            out_seq.append(SeqIO.SeqRecord(Seq(seq_str), id=key))
        
        out_fasta_name = os.path.join(out_dir, '_'.join(sample_list) + "_" + geneid + ".fa")

        with open(out_fasta_name, 'w') as f:
            pass  # Clear the file
        
        # Write results to file
        with open(out_fasta_name, "w") as output_handle:
            SeqIO.write(out_seq, output_handle, "fasta")

