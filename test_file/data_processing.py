import scanpy as sc

def process_data(file_path):
    adata = sc.read_h5ad(file_path)
    
    return adata