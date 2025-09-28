import pandas as pd
import os

def get_id(protein_id):
    return protein_id.split("-")[0]

def count_n_domain_gene(domain, gene_lis, interpro_data, message = True):
    """
    domain: string signature_accession 

    gene_lis: list, 
    """
    n_domain_gene = 0
    for gene in gene_lis:
        df = show_interpro_feature_lis([gene], interpro_data)
        if df["signature_accession"].isin([domain]).any():
            #print(df[df["signature_accession"] == "TRANSMEMBRANE"])
            if message:
                print(f"{gene} contains {domain} signature")
            n_domain_gene += 1
        else:
            if message:
                print(f"{gene} does not contain {domain} signature")
    if message: 
        print(f'''\
        Number of genes with {domain}: {n_domain_gene}
        Number of total genes: {len(gene_lis)}''')
    return n_domain_gene
    
def count_n_domain_iprid_gene(domain, iprid, KaKs_data, high_kaks_gene_col, interpro_data):
    genes = KaKs_data[KaKs_data["interpro_accession"] == iprid][high_kaks_gene_col].iloc[0]
    gene_lis = list(set(genes.split(",")))
    n_gene = len(gene_lis)
    n_domain_gene = count_n_domain_gene(domain, gene_lis, interpro_data, message = False)
    return (n_gene, n_domain_gene)
def show_interpro_feature_lis(geneid_lis, interpro_data):
    pd.set_option('display.max_rows', 500)
    return interpro_data[interpro_data["geneid"].isin(geneid_lis)]

def count_n_ipr_acc_iprid_gene(ipr_acc, iprid_gene, KaKs_data, high_kaks_gene_col, interpro_data):
    genes = KaKs_data[KaKs_data["interpro_accession"] == iprid_gene][high_kaks_gene_col].iloc[0]
    gene_lis = list(set(genes.split(",")))
    n_gene = len(gene_lis)
    n_domain_gene = 0
    for gene in gene_lis:
        df = show_interpro_feature_lis([gene], interpro_data)
        if df["interpro_accession"].isin([ipr_acc]).any():
            #print(df[df["signature_accession"] == "TRANSMEMBRANE"])
            n_domain_gene += 1
    return (n_gene, n_domain_gene)
def count_genes(vec):
    vec = str(vec)
    return len(set(vec.split(",")))
def get_gene_set(vec):
    vec = str(vec)
    gene_lis = list(set(vec.split(",")))
    gene_lis.sort()
    return ",".join(gene_lis)


def get_pair_w_most_nonsyn(geneid, result_folder, samples = ["C1G11", "C1C5v2", "C3D7v2"]):
    #geneid = "FUN_000296"
    file_name = geneid + "_NonSyn_count.csv"
    file_path = os.path.join(result_folder, file_name)
    df = pd.read_csv(file_path)
    df_subset = df.loc[1:,samples].reset_index(drop = True)
    result_dict = {}
    for sampleid in range(len(samples)):
        for sampleid2 in range(sampleid,len(samples)):
            sample1 = samples[sampleid]
            #print(sample1)
            sample2 = samples[sampleid2]
            #print(sample2)
            key = sample1 + "_" + sample2
            result_dict[key] = df_subset.loc[sampleid, sample2]
    if max(result_dict.values()) == 0:
        return "No_NonSyn_substitutions"
    else:
        return max(result_dict, key = result_dict.get)