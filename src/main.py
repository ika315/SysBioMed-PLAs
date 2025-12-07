import os
import utils

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
gene_csv = os.path.join(base_dir, "data", "updated_gene_list.csv")

if __name__ == "__main__":
    genes = utils.read_gene_list(gene_csv)
    print(genes)

