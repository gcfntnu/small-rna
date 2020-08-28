#!/usr/bin env python

import sys
import argparse
import warnings
warnings.filterwarnings('ignore', message='numpy.dtype size changed')

import pandas as pd
import numpy as np
import scanpy as sc


def argparser():
    parser = argparse.ArgumentParser(description='Dimred plot from anndata')
    parser.add_argument('-i', '--input',
                        help='Input filename. Anndata file (.h5ad)')
    parser.add_argument('-o ', '--output', default='umap_mqc.png',
                        help='Output filename. Will default to pca_mqc.png, Optional')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = argparser()
    
    adata = sc.read(args.input)

    sc.pp.filter_cells(adata, min_genes=10)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.1, max_mean=20, min_disp=0.2)

    n_comps = min(20, min(adata.shape)-1)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
    #sc.pp.neighbors(adata, n_neighbors=min(10, n_comps), n_pcs=n_comps)
    #tl.umap(adata)
    #sc.pl.pca(adata, color=['Sample_Group'])
    fig = sc.pl.pca(adata, color=['Sample_Group'], return_fig=True)
    fig.savefig(args.output)
