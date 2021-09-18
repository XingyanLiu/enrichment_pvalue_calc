# -*- coding: UTF-8 -*-
"""
@Author: Xingyan Liu
@CreateDate: 2021-09-12
@File: run_analysis.py
@Project: pvalue_calc
"""
import os
import sys
from pathlib import Path
import logging
import pandas as pd
import numpy as np
sys.path.append('.')
from hyperg_test_pvalue import load_namelist, multi_sets_pvalues, enrichment_pvalue


def _top_n_markers(df: pd.DataFrame, ntop: int = 10):
    return pd.unique(df.head(ntop).values.flatten())


def _load_gene_universe() -> set:
    datadir = Path(__file__).parent
    return set(load_namelist(datadir / "genes-dog.tsv"))


def _load_cluster_deg_dict() -> dict:
    datadir = Path(__file__).parent
    fn = datadir / 'DEGs-dog/deg_intersects-p0.001.tsv'
    degs_each_cl = pd.read_csv(fn, sep='\t', index_col=0).iloc[:, 0]
    return degs_each_cl.map(lambda x: x.split(',')).to_dict()


def __load_pooled_cluster_deg() -> set:
    deg = set()
    deg_dict = _load_cluster_deg_dict()
    for cl, s in deg_dict.items():
        deg.update(s)
    return deg


def _load_psg_dict() -> dict:
    datadir = Path(__file__).parent / 'ref_sets'
    psg_tags = ["plos2016", "moreThan2", "moreThan3", "0"]
    return {t: set(load_namelist(datadir / f"psg-{t}.txt")) for t in psg_tags}


def _load_topk_deg_dict(de_test='MAST') -> dict:
    datadir = Path(__file__).parent
    genes_mk_dict = {}
    for _n in [10, 20, 30, 50]:
        file_markers = datadir / f"DEGs-dog/top-{_n}-cluster({de_test}).csv"
        genes_mk = load_namelist(file_markers)
        genes_mk_dict[f"top-{_n}"] = set(genes_mk)
    return genes_mk_dict


def main_0():
    datadir = Path(__file__).parent
    resdir = datadir / 'pvalues'
    if not os.path.exists(resdir):
        os.mkdir(resdir)
    genes_all = _load_gene_universe()
    # PSGs
    psg_dict = _load_psg_dict()
    # DEGs
    de_test = "MAST"
    for de_test in ['MAST', 't', 'wilcox']:
        genes_mk_dict = _load_topk_deg_dict(de_test)

        resdf = multi_sets_pvalues(psg_dict, genes_mk_dict, genes_all)
        print(resdf)
        resdf.to_csv(resdir / f'pvalues-DEby_{de_test}.csv', index_label="ntop_markers")


def main_each_cluster():
    datadir = Path(__file__).parent
    resdir = datadir / 'pvalues'
    if not os.path.exists(resdir):
        os.mkdir(resdir)
    # universe genes
    genes_all = _load_gene_universe()
    # PSGs
    psg_dict = _load_psg_dict()
    # DEGs
    degs_each_cl = _load_cluster_deg_dict()

    resdf = multi_sets_pvalues(psg_dict, degs_each_cl, genes_all)
    print(resdf)
    resdf.to_csv(resdir / f'pvalues-perCluster-IntersectDEGs.csv', index_label="cluster")


def main_pooled_degs():
    datadir = Path(__file__).parent
    resdir = datadir / 'pvalues'
    if not os.path.exists(resdir):
        os.mkdir(resdir)
    # universe genes
    genes_all = _load_gene_universe()
    # PSGs
    psg_dict = _load_psg_dict()
    # DEGs
    degs_each_cl = __load_pooled_cluster_deg()
    resdf = multi_sets_pvalues(
        psg_dict, {'pooled DEGs': degs_each_cl},
        genes_all)

    print(resdf)
    resdf.to_csv(resdir / f'pvalues-all-IntersectDEGs.csv', )


if __name__ == "__main__":
    # datadir = Path("/Users/xingyan/Downloads/temp/pvalue_calc/")
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(filename)s-%(lineno)d-%(funcName)s(): '
               '%(levelname)s\n %(message)s')
    # main_0()
    # main_each_cluster()
    main_pooled_degs()
