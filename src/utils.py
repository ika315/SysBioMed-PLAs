import pandas as pd

""" possible plan:
- take curated gene list
- filter for genes detectable in selected data set?
- compute PLA score
- rank genes based on correlation with PLA score?
- maybe correlation between genes expression and PLA score?
- choose genes that are above a certain threshold 
- outcome = a refined set of biologically meaningful candidate PLA-associated genes?
- update gene set
"""

def read_gene_list(gene_csv):
    use_cols = ["geneName"]
    gene_list = pd.read_csv(gene_csv, usecols=use_cols)
    genes = set(gene_list["geneName"])
    print(f"{len(genes)} genes were extracted from file.")
    return genes

def extend_gene_set(gene_set):
    pass