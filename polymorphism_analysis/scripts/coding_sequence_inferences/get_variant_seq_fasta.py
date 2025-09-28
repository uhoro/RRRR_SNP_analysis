import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
import os
import argparse
from os import listdir
from os.path import isfile, join


def get_sample_sequence_fast(filepath, sampleid):
    with open(filepath, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if sampleid in record.id:
                return SeqIO.SeqRecord(Seq(str(record.seq)), id=sampleid)
    return None

def get_variant_seq_fa(variant_seq_folder, 
                    sample_list, sample_dir_list, 
                    test = False, include_ref = True):

    os.makedirs(variant_seq_folder, exist_ok=True)

    for j in range(len(sample_list)):
        print(f"sample id: {sample_list[j]}")
        sample_files = [f for f in listdir(sample_dir_list[j]) if isfile(join(sample_dir_list[j], f))]
        print(f"Number of seqs in {sample_dir_list[j]}: {len(sample_files)}")

    sample1_files = [f for f in listdir(sample_dir_list[0]) if isfile(join(sample_dir_list[0], f))]

    # for test, do only for the first 10 geneid_vec
    if test:
        last = 10
    else:
        last = len(sample1_files)
    
    i = 0
    for file in sample1_files:
        if i >= last:
            break
        i += 1
        if i % 100 == 0:
            print(f"{i} th gene")
            print(f"start processing: {gene_name}")

        out_sequences = []
        if include_ref:
            # fasta_sequences1 = list(SeqIO.parse(open(sample_dir_list[0] + "/" + file),'fasta'))
            # ref_rec = [sequence for sequence in fasta_sequences1 if "REF" in sequence.id]
            # ref_rec_new = SeqIO.SeqRecord(
            #     Seq(ref_rec[0].seq),
            #     id = "REF"
            # )
            # out_sequences.append(ref_rec_new)
            ref_rec = get_sample_sequence_fast(sample_dir_list[0] + "/" + file, "REF")
            out_sequences.append(ref_rec)


        for sampleid in range(0, len(sample_list)):
            # sample_files = [f for f in listdir(sample_dir_list[sampleid]) if isfile(join(sample_dir_list[sampleid], f))]

            # if file not in sample_files:
            #     print(f"{file} not found in {sample_dir_list[sampleid]} folder")
            #     continue

            # fasta_sequences = list(SeqIO.parse(open(sample_dir_list[sampleid] + "/" + file),'fasta'))
            # sample_rec = [sequence for sequence in fasta_sequences if sample_list[sampleid] in sequence.id]
            # sample_rec_new = SeqIO.SeqRecord(
            #     Seq(sample_rec[0].seq),
            #     id = sample_list[sampleid]
            # )
            # out_sequences.append(sample_rec_new)
            sample_rec = get_sample_sequence_fast(sample_dir_list[sampleid] + "/" + file, sample_list[sampleid])
            out_sequences.append(sample_rec)

        gene_id = file.split(".")[0]
        gene_name = gene_id.split("-")[0]

        ref = ""
        if include_ref:
            ref = "REF_"
        out_fasta_name = os.path.join(variant_seq_folder, ref + '_'.join(sample_list) + "_" + gene_name + ".fa")

        with open(out_fasta_name, 'w') as f:
            pass  # Clear the file
        
        # Write results to file
        with open(out_fasta_name, "w") as output_handle:
            SeqIO.write(out_sequences, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get in fasta format variant sequences of entered sample list from vcf2fasta result"
    )
    parser.add_argument("--out_dir", required=True, type=str, 
                        help = "Enter the directory where variant sequences should be stored")
    parser.add_argument("--sample_list", nargs='+', default = ["C1G11", "C1C5v2", "C3D7v2"],
                        help = "Enter sample list. Default is C1G11 C1C5v2 C3D7v2")
    parser.add_argument("--sample_dir_list", nargs='+', default = ["out/set250615_C1G11_CDS", 
                                                                   "out/set250615_C1C5v2_CDS",
                                                                   "out/set250615_C3D7v2_CDS"],
                        help = "Enter sample directory list. Default is out/set250615_C1G11_CDS" \
                                                                   "out/set250615_C1C5v2_CDS" \
                                                                   "out/set250615_C3D7v2_CDS")
    parser.add_argument("--test", required=False, type=bool, default = False, 
                        help = "Enable test mode with True")
    parser.add_argument("--include_ref", required=False, type=bool, default = True, 
                        help = "Remove reference sequence by setting False")
    args = parser.parse_args()

    variant_seq_folder = args.out_dir
    sample_list = args.sample_list
    sample_dir_list = args.sample_dir_list
    test = args.test
    include_ref = args.include_ref
    #print("sample_list")
    #print(sample_list)
    #print(type(sample_list))


    get_variant_seq_fa(variant_seq_folder,
                   sample_list, sample_dir_list, 
                   test = test, include_ref = include_ref)












