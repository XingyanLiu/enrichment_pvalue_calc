# -*- coding: UTF-8 -*-
"""
The hypergeometric distribution models drawing objects from a bin.
M is the total number of objects, n is total number of Type I objects.
The random variate represents the number of Type I objects in N drawn without
replacement from the total population.

Apply to gene enrichment analysis

* M = Gene_universe_size
* n = Ontology_gene_set_size
* N = Selected_gene_set_size
* k = Intersected_gene_set
"""
import os
import random
from pathlib import Path
from typing import Sequence, Union, Dict, List, Set, Any
import pandas as pd
from scipy.stats import hypergeom


def load_namelist(filename):
    with open(filename) as f:
        lines = filter(lambda x: len(x.strip()) >= 1, f.readlines())
        return list(map(lambda x: x.strip().strip('"'), lines))


def intersect(set1, set2):
    return set(set1).intersection(set2)


def hypergeom_pvalue(
    k: int or Sequence[int],  # intersection
    N: int or Sequence[int],  # selected gene set size(s), e.g. marker genes
    n: int,  # PSG size
    M: int,  # Gene universe size
):
    if isinstance(k, int):
        return hypergeom.sf(k, M, n, N)
    foo_pval = lambda k, N: hypergeom.sf(k, M, n, N)
    return list(map(foo_pval, k, N))


def pvalue_random(
        n_random: int,
        universe: Sequence or Set,
        ref_set: Sequence or Set,
        n_repeats: int = 1000,
):
    """ Compute the background p-value

    :param n_random:
    :param universe:
    :param ref_set:
        You should make sure that ALL the elements in ref_set exists in the universe-set
    :param n_repeats:
        Number of repeats
    :return:
    """
    pvalues = 0
    for i in range(n_repeats):
        random_set = random.sample(universe, n_random)
        k = len(intersect(random_set, ref_set))
        M, n = len(universe), len(ref_set)
        pvalues += hypergeom_pvalue(k, n_random, n, M)
    return pvalues / n_repeats


def enrichment_pvalue(
        ref_set: Sequence or Set,
        que_set: Sequence or Set,
        universe: Union[List, Set],
        check_valid: bool = True,
):
    """ Compute the enrichment p-value based on hypergeometric distribution models

    :param ref_set:
        referent set of elements (that is known to be of interest)
    :param que_set:
        query set of elements, e.g. DEGs
    :param universe:
        the background elements
    :param check_valid:
        whether to filter out the elements that are not in the universe
    :return:

    for a list or dict of que-sets:
    >>> pvals = [enrichment_pvalue(ref_set, s) for s in que_sets]
    >>> pvals = {k: enrichment_pvalue(ref_set, s) for k, s in que_sets.items()}

    """
    universe = set(universe)
    if check_valid:
        ref_set = set(x for x in ref_set if x in universe)
        que_set = set(x for x in que_set if x in universe)
    k = len(intersect(que_set, ref_set))
    M, n, N = len(universe), len(ref_set), len(que_set)
    pval = hypergeom_pvalue(k, N, n, M)
    # print(f"{pval:.3f}")
    return pval


def multi_sets_pvalues(
        ref_sets: Union[Sequence[Set], Dict[Any, Set]],  # PSGs
        que_sets: Union[Sequence[Set], Dict[Any, Set]],  # markers
        universe: Union[List, Set],
        n_rep_random: int = 1000,
) -> pd.DataFrame:
    """ Packed function for multiple reference and query sets

    :param ref_sets: list or dict of sets
    :param que_sets: list or dict of sets
    :param universe: background set or list
    :param n_rep_random:
        number of random repeats (as background comparison)
    :return: pd.DataFrame
        each column for a reference and each row for a query set (or random).
    example:
    >>> resdf = multi_sets_pvalues(psg_dict, degs_each_cl, genes_all)
    >>> print(resdf)
    >>> resdf.to_csv(resdir / f'pvalues-perCluster.csv', index_label="cluster")

    """
    if isinstance(ref_sets, Sequence):
        ref_tags = [f'ref_{i}' for i in range(len(ref_sets))]
        ref_sets = dict(zip(ref_tags, ref_sets))
    if isinstance(que_sets, Sequence):
        que_tags = [f'que_{i}' for i in range(len(que_sets))]
        que_sets = dict(zip(que_tags, que_sets))
    ref_sets = {k: intersect(s, universe) for k, s in ref_sets.items()}
    que_sets = {k: intersect(s, universe) for k, s in que_sets.items()}

    record = {}
    for ref_tag, ref_set in ref_sets.items():
        record[ref_tag] = {}
        for q_tag, qset in que_sets.items():
            record[ref_tag][q_tag] = enrichment_pvalue(
                ref_set, qset, universe, check_valid=False)
            # random as background
            record[ref_tag][f"{q_tag}(random)"] = pvalue_random(
                len(qset), universe, ref_set, n_repeats=n_rep_random)

    resdf = pd.DataFrame(record)
    return resdf


def __test__():
    datadir = Path(__file__).parent / 'test_data'
    universe = load_namelist(datadir / 'universeSet.txt')
    ref_set = load_namelist(datadir / 'refSet.txt')
    que_set = load_namelist(datadir / 'querySet.txt')

    # test 1
    res = enrichment_pvalue(ref_set, que_set, universe)
    logging.info(res)

    # test 2
    res = multi_sets_pvalues(
        [ref_set] * 2, [que_set] * 2, universe, n_rep_random=1000)
    logging.info(res)


if __name__ == "__main__":
    import logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(filename)s-%(lineno)d-%(funcName)s(): '
               '%(levelname)s\n %(message)s')

    __test__()


