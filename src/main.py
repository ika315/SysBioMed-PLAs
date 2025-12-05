import os
import pandas as pd

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
gene_csv = os.path.join(base_dir, "data", "platelet-genes.csv")

def read_gene_list(gene_csv):
    use_cols = ["geneName"]
    gene_list = pd.read_csv(gene_csv, usecols=use_cols)
    genes = set(gene_list["geneName"])
    print(len(genes))
    #print(gene_list.head(n=10))

if __name__ == "__main__":
    read_gene_list(gene_csv)

