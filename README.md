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

