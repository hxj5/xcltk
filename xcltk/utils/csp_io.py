# csp_io.py - utils for processing cellsnp-lite input and output.


import anndata as ad
import gzip
import numpy as np
import os
import pandas as pd
import pysam
import scipy as sp
from scipy import io
from scipy import sparse



def load_data(data_dir, is_genotype = False, is_gzip = True):
    vcf_suffix = ".gz" if is_gzip else ""
    base_vcf_comment, base_vcf = load_vcf(
        os.path.join(data_dir, "cellSNP.base.vcf" + vcf_suffix))
    base_vcf = base_vcf.rename(columns = {
        'CHROM': 'chrom',
        'POS': 'pos',
        'REF': 'ref',
        'ALT': 'alt'
    })
    
    cell_vcf, cell_vcf_comment = None, None
    if is_genotype:
        cell_vcf_comment, cell_vcf = load_vcf(
            os.path.join(data_dir, "cellSNP.cells.vcf" + vcf_suffix))
        cell_vcf = cell_vcf.rename(columns = {
            'CHROM': 'chrom',
            'POS': 'pos',
            'REF': 'ref',
            'ALT': 'alt'
        })

    cells = load_cells(os.path.join(data_dir, "cellSNP.samples.tsv"))
    AD_mtx = load_matrix(os.path.join(data_dir, "cellSNP.tag.AD.mtx"))
    DP_mtx = load_matrix(os.path.join(data_dir, "cellSNP.tag.DP.mtx"))
    OTH_mtx = load_matrix(os.path.join(data_dir, "cellSNP.tag.OTH.mtx"))
    
    adata = ad.AnnData(
        X = None, 
        obs = base_vcf, 
        var = cells
    )
    
    adata.uns["is_gzip"] = is_gzip
    adata.uns["is_genotype"] = is_genotype
    adata.uns["base_vcf_comment"] = base_vcf_comment
    if cell_vcf is None:
        adata.uns["cell_vcf_comment"] = None
    else:
        adata.obsm["cell_vcf"] = cell_vcf
        adata.uns["cell_vcf_comment"] = cell_vcf_comment
    
    adata.layers["AD"] = AD_mtx
    adata.layers["DP"] = DP_mtx
    adata.layers["OTH"] = OTH_mtx
    
    adata = adata.transpose()        # to yield cell x feature adata.
    return(adata)


def save_data(adata, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    is_gzip = adata.uns["is_gzip"]
    is_genotype = adata.uns["is_genotype"]
    
    adata = adata.transpose()        # to yield feature x cell adata.
    
    save_cells(adata.var, os.path.join(out_dir, "cellSNP.samples.tsv"))
    vcf_suffix = ".gz" if is_gzip else ""
    save_vcf(
        adata.obs, 
        comment = adata.uns["base_vcf_comment"], 
        fn = os.path.join(out_dir, "cellSNP.base.vcf" + vcf_suffix),
        is_gzip = is_gzip
    )
    if is_genotype:
        save_vcf(
            adata.obsm["cell_vcf"],
            comment = adata.uns["cell_vcf_comment"],
            fn = os.path.join(out_dir, "cellSNP.cells.vcf" + vcf_suffix),
            is_gzip = is_gzip
        )

    save_matrix(adata.layers["AD"], os.path.join(out_dir, "cellSNP.tag.AD.mtx"))
    save_matrix(adata.layers["DP"], os.path.join(out_dir, "cellSNP.tag.DP.mtx"))
    save_matrix(adata.layers["OTH"], os.path.join(out_dir, "cellSNP.tag.OTH.mtx"))
    
    

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



def load_cells(fn):
    df = pd.read_csv(fn, header = None)
    df.columns = ["cell"]
    return(df)


def save_cells(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)



def load_vcf(fn):
    # load comment
    fp = None
    if fn.endswith(".gz") or fn.endswith(".GZ"):
        fp = gzip.open(fn, "rt")
    else:
        fp = open(fn, "r")
        
    comment = ""
    pre_line = None
    for line in fp:
        if not line or line[0] != "#":
            break
        pre_line = line
        comment += line
    
    fp.close()
    
    if not pre_line:
        raise IOError()
    assert len(pre_line) > 6
    assert pre_line[:6] == "#CHROM"
    
    # load content
    content = pd.read_csv(fn, sep = "\t", header = None, comment = "#")
    content.columns = pre_line.strip()[1:].split("\t")
    content['CHROM'] = content['CHROM'].astype(str)
    
    return((comment, content))


def save_vcf(df, comment, fn, is_gzip = True):
    df_file = fn + ".df"
    df.to_csv(df_file, sep = "\t", header = False, index = False)
    
    fp = pysam.BGZFile(fn, "w") if is_gzip else open(fn, "w")
    fp.write(comment.encode())
    with open(df_file, "r") as df_fp:
        for line in df_fp:
            fp.write(line.encode())
    fp.close()
    
    os.remove(df_file)
