# io.py - xcltk RDR module input and output.

import anndata as ad
import gzip
import numpy as np
import os
import pandas as pd
import scipy as sp
from scipy import io
from scipy import sparse


    
def load_data(data_dir):
    features = load_features(os.path.join(data_dir, "features.tsv"))
    cells = load_cells(os.path.join(data_dir, "barcodes.tsv"))
    mtx = load_matrix(os.path.join(data_dir, "matrix.mtx"))
    
    adata = ad.AnnData(
        X = mtx, 
        obs = features,
        var = cells
    )
    
    adata = adata.transpose()        # to yield cell x feature adata.
    return(adata)


def save_data(adata, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    save_cells(adata.obs, 
        fn = os.path.join(out_dir, "barcodes.tsv")) 
    save_features(adata.var, 
        fn = os.path.join(out_dir, "features.tsv"))
    save_matrix(adata.X, os.path.join(out_dir, "matrix.mtx"))


    
def load_matrix(fn):
    mtx = None
    try:
        mtx = sp.io.mmread(fn)
    except:
        mtx = io.mmread(fn)
    mtx = mtx.toarray()    # convert from sparse matrix to ndarray to support slicing.
    return(mtx)


def save_matrix(mtx, fn):
    mtx = sparse.csr_matrix(mtx)   # convert from ndarray to sparse matrix to be fully compatible with .mtx format.
    io.mmwrite(fn, mtx)

    

def load_features(fn):
    df = pd.read_csv(fn, header = None, sep = "\t", dtype = {0: str})
    df.columns = ["chrom", "start", "end", "feature"]
    return(df)


def save_features(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


    
def load_cells(fn):
    df = pd.read_csv(fn, header = None)
    df.columns = ["cell"]
    return(df)


def save_cells(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)
