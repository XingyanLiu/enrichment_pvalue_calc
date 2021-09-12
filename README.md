# Compute the enrichment p-value based on hypergeometric distribution models

The hypergeometric distribution models drawing objects from a bin.
M is the total number of objects, n is total number of Type I objects.
The random variate represents the number of Type I objects in N drawn without
replacement from the total population.

When applied to gene-enrichment analysis:

* M = Gene_universe_size
* n = Ontology_gene_set_size
* N = Selected_gene_set_size
* k = Intersected_gene_set

Test example code:

```python
from pathlib import Path
from hyperg_test_pvalue import load_namelist, enrichment_pvalue, multi_sets_pvalues

datadir = Path('./test_data')
universe = load_namelist(datadir / 'universeSet.txt')
ref_set = load_namelist(datadir / 'refSet.txt')
que_set = load_namelist(datadir / 'querySet.txt')

# test 1
print(" TEST 1 ".center(60, '='))
res = enrichment_pvalue(ref_set, que_set, universe)
print(res)

# test 2
print(" TEST 2 ".center(60, '='))
res = multi_sets_pvalues(
    [ref_set] * 2, [que_set] * 2, universe, n_rep_random=1000)
print(res)
```