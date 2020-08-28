#!/usr/bin env python

import sys
import argparse
import warnings
warnings.filterwarnings('ignore', message='numpy.dtype size changed')

import pandas as pd
import numpy as np
import anndata as ad


def argparser():
    parser = argparse.ArgumentParser(description='Create anndata file from tab counts')
    parser.add_argument('-i', '--input',
                        help='Input filename. Tab separated mirna counts.')
    parser.add_argument('--sample-info', dest='samples',
                        help='Optional sample sheet. Will subset expr table if needed')
    parser.add_argument('--feature-info', dest='features',
                        help='Optional sample sheet. Will subset expr table if needed')
    parser.add_argument('-o ', '--output', default='adata.h5ad',
                        help='Output filename. Will default to pca_mqc.png, Optional [*.h5ad]')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = argparser()
    
    X = pd.read_csv(args.input, sep="\t", index_col=0).T
    X.columns = X.columns.astype(str)

    S = F = None
    if args.samples is not None:
        S = pd.read_csv(args.samples, sep="\t", index_col=0)
        S.index = S.index.astype(str)
        if not X.index.isin(S.index).all():
            print(S.head())
            print(X.columns)
            raise ValueError("missing samples in sample info!")
        S = S.loc[X.index, :]
        
    if args.features is not None:
        F = pd.read_csv(args.features, sep="\t", index_col=0)
        if not X.columns.isin(F.index).all():
            warnings.warn("missing annotations in feature info!")
            F = F.reindex(X.columns)
            
        F = F.loc[X.columns,:]

    adata = ad.AnnData(X=X, obs=S, var=F)

    adata.write(filename=args.output)
    
