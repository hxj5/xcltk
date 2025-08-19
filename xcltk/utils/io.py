# io.py - xcltk input and output.

import anndata as ad
import gzip
import numpy as np
import os
import pandas as pd
import scipy as sp
from scipy import io
from scipy import sparse



def baf_load_data(data_dir):
    regions = load_regions(os.path.join(data_dir, "xcltk.region.tsv"))
    samples = load_samples(os.path.join(data_dir, "xcltk.samples.tsv"))
    AD_mtx = load_matrix(os.path.join(data_dir, "xcltk.AD.mtx"))
    DP_mtx = load_matrix(os.path.join(data_dir, "xcltk.DP.mtx"))
    OTH_mtx = load_matrix(os.path.join(data_dir, "xcltk.OTH.mtx"))
    
    adata = ad.AnnData(
        X = None, 
        obs = regions,
        var = samples)
    
    adata.layers["AD"] = AD_mtx
    adata.layers["DP"] = DP_mtx
    adata.layers["OTH"] = OTH_mtx
    
    return(adata)


def baf_save_data(adata, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    save_regions(adata.obs, 
        fn = os.path.join(out_dir, "xcltk.region.tsv"))
    
    save_samples(adata.var, 
        fn = os.path.join(out_dir, "xcltk.samples.tsv"))
    
    save_matrix(adata.layers["AD"], os.path.join(out_dir, "xcltk.AD.mtx"))
    save_matrix(adata.layers["DP"], os.path.join(out_dir, "xcltk.DP.mtx"))
    save_matrix(adata.layers["OTH"], os.path.join(out_dir, "xcltk.OTH.mtx"))


    
def rdr_load_data(data_dir):
    regions = load_regions(os.path.join(data_dir, "features.tsv"))
    samples = load_samples(os.path.join(data_dir, "barcodes.tsv"))
    mtx = load_matrix(os.path.join(data_dir, "matrix.mtx"))
    
    adata = ad.AnnData(
        X = mtx, 
        obs = regions,
        var = samples)
    
    return(adata)


def rdr_save_data(adata, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    save_regions(adata.obs, 
        fn = os.path.join(out_dir, "features.tsv"))
    
    save_samples(adata.var, 
        fn = os.path.join(out_dir, "barcodes.tsv"))
    
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

    

def load_regions(fn):
    df = pd.read_csv(fn, header = None, sep = "\t", dtype = {0: str})
    df.columns = ["chrom", "start", "end", "feature"]
    return(df)


def save_regions(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)


    
def load_samples(fn):
    df = pd.read_csv(fn, header = None)
    df.columns = ["cell"]
    return(df)


def save_samples(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)
