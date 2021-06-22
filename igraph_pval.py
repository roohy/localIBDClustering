import numpy as np
from scipy.stats import binom_test


def pval_from_igraph_partitions(parition,pheno_preval,pheno_array):
    partition_subgraphs = parition.subgraphs()
    pvals = np.zeros(len(pheno_preval))
    for subgraph in partition_subgraphs:
        result_dict = {}
        for vertex in subgraph.vs:
            # phenos = np.where(pheno_array[:,vertex['name']])[0]
            phenos = pheno_array[int(vertex['name'])]
            if len(phenos) != 0:
                for pheno in phenos:
                    if pheno in result_dict:
                        result_dict[pheno] += 1
                    else:
                        result_dict[pheno] = 1
        cls_size = len(subgraph.vs)
        for key in result_dict:
            pval = binom_test(result_dict[key],cls_size,pheno_preval[key])
            pval = np.abs(np.log10(pval))
            if pvals[key] < pval:
                pvals[key] = pval
    return pvals