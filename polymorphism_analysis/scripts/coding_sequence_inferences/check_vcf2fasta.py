import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import argparse


def compare_sequences(seq, ref):
    if len(seq) == len(ref):
        if seq == ref:
            return "identical"
        else:
            return "same length, different sequences"
    else:
        if "-" in seq:
            seq_wo_gaps = seq.replace("-", "")
            if seq_wo_gaps == ref:
                return "identical when gaps are removed"
            else:
                return "different sequences after gaps removed"
        else:
            return "different lengths and no gaps"



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get gap information in the fasta files in given directory"
    )
    parser.add_argument("--out_dir", required=True, type=str, 
                        help = "Enter the directory where the failed cds gene information data is stored")
    parser.add_argument("--out_fname", required=False, type=str, default = "failed_cds.txt",
                        help = "Enter the file name storing the failed cds gene information")
    parser.add_argument("--msa_data_dir", required=True, type=str, 
                        help = "Enter sample id for which gap information is retrieved")
    parser.add_argument("--ref_cds_file", required = True, type=str,
                        help = "Enter the directory of original fasta file")

    args = parser.parse_args()

    out_dir = args.out_dir
    out_fname = args.out_fname
    msa_data_dir = args.msa_data_dir
    ref_cds_file = args.ref_cds_file

    os.makedirs(out_dir, exist_ok=True)

    ref_sequences = {}
    for record in SeqIO.parse(ref_cds_file , "fasta"):
        record_id = record.id.split("-")[0]
        ref_sequences[record_id] = str(record.seq)
    gene_list = list(ref_sequences.keys())

    result_array = []
    i = 0
    for geneid in gene_list:
        i += 1
        if i % 1000 == 0:
            print(f"Processing gene {i}/{len(gene_list)}: {geneid}")
        msa_fname = "REF_C1G11_C1C5v2_C3D7v2_" + geneid + ".fa"
        msa_data_full_path = os.path.join(msa_data_dir, msa_fname)
        file = Path(msa_data_full_path)
        if file.is_file():
            msa_seqs = {}
            for record in SeqIO.parse(msa_data_full_path , "fasta"):
                msa_seqs[record.id] = str(record.seq)
            result_array.append(compare_sequences(msa_seqs["REF"], ref_sequences[geneid]))
        else:
            result_array.append("msa file not found")

    ref_validation_df = pd.DataFrame({
                            "geneid": gene_list,
                            "validation_result": result_array
                        })
 
    failed_cds = ref_validation_df[ref_validation_df["validation_result"] == "different lengths and no gaps"]["geneid"].tolist() +  \
                    ref_validation_df[ref_validation_df["validation_result"] == "different sequences after gaps removed"]["geneid"].tolist() + \
                    ref_validation_df[ref_validation_df["validation_result"] == "msa file not found"]["geneid"].tolist()
    
    with open(os.path.join(out_dir, out_fname), "w") as f:
        for gene in failed_cds:
            line = gene + "\n"
            f.write(line)

