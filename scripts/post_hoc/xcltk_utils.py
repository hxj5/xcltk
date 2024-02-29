# xcltk_utils.py - utils for processing xcltk input and output.

import anndata as ad
import gzip
import numpy as np
import os
import pandas as pd
import scipy as sp
from scipy import io
from scipy import sparse


def xtk_baf_load_data(data_dir):
    regions = xtk_load_regions(os.path.join(data_dir, "xcltk.regions.tsv"))
    samples = xtk_load_samples(os.path.join(data_dir, "xcltk.samples.tsv"))
    AD_mtx = xtk_load_matrix(os.path.join(data_dir, "xcltk.AD.mtx"))
    DP_mtx = xtk_load_matrix(os.path.join(data_dir, "xcltk.DP.mtx"))
    OTH_mtx = xtk_load_matrix(os.path.join(data_dir, "xcltk.OTH.mtx"))
    
    adata = ad.AnnData(
        X = AD_mtx, 
        obs = regions,
        var = samples)
        
    adata.layers["DP"] = DP_mtx
    adata.layers["OTH"] = OTH_mtx
    
    return(adata)


def xtk_baf_save_data(adata, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    xtk_save_regions(adata.obs, 
        fn = os.path.join(out_dir, "xcltk.regions.tsv"))
    
    xtk_save_samples(adata.var, 
        fn = os.path.join(out_dir, "xcltk.samples.tsv"))
    
    xtk_save_matrix(adata.X, os.path.join(out_dir, "xcltk.AD.mtx"))
    xtk_save_matrix(adata.layers["DP"], os.path.join(out_dir, "xcltk.DP.mtx"))
    xtk_save_matrix(adata.layers["OTH"], os.path.join(out_dir, "xcltk.OTH.mtx"))


def xtk_rdr_load_data(data_dir):
    regions = xtk_load_regions(os.path.join(data_dir, "features.tsv"))
    samples = xtk_load_samples(os.path.join(data_dir, "barcodes.tsv"))
    mtx = xtk_load_matrix(os.path.join(data_dir, "matrix.mtx"))
    
    adata = ad.AnnData(
        X = mtx, 
        obs = regions,
        var = samples)
    
    return(adata)


def xtk_rdr_save_data(adata, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    xtk_save_regions(adata.obs, 
        fn = os.path.join(out_dir, "features.tsv"))
    
    xtk_save_samples(adata.var, 
        fn = os.path.join(out_dir, "barcodes.tsv"))
    
    xtk_save_matrix(adata.X, os.path.join(out_dir, "matrix.mtx"))


def xtk_load_matrix(fn):
    mtx = None
    try:
        mtx = sp.io.mmread(fn)
    except:
        mtx = io.mmread(fn)
    mtx = mtx.toarray()    # convert from sparse matrix to ndarray to support slicing.
    return(mtx)


def xtk_load_regions(fn):
    df = pd.read_csv(fn, header = None, sep = "\t")
    df.columns = ["chrom", "start", "end", "reg_id"]
    return(df)


def xtk_load_samples(fn):
    df = pd.read_csv(fn, header = None)
    df.columns = ["cell"]
    return(df)


def xtk_save_matrix(mtx, fn):
    mtx = sparse.csr_matrix(mtx)   # convert from ndarray to sparse matrix to be fully compatible with .mtx format.
    io.mmwrite(fn, mtx)


def xtk_save_regions(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


def xtk_save_samples(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


#if __name__ == "__main__":
#    pass