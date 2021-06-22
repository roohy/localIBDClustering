import numpy as np
from scipy.stats import binom_test


'''def test_igraph_partitions(partition,pheno_preval,pheno_array):
    partition_subgraphs = partition.subgraphs()
    pvals = np.zeros(len(pheno_preval))
    for subgraph in partition_subgraphs:
        result_dict = {}
        for vertext in subgraph.vs:
            if pheno_array[vertext['name']] != -1:
                if pheno_array[vertext['name']] in result_dict:
                    result_dict[pheno_array[vertext['name']]] += 1
                else:
                    result_dict[pheno_array[vertext['name']]]  = 1
        cls_size = len(subgraph.vs)
        for key in result_dict:
            pval = binom_test(result_dict[key],cls_size,pheno_preval[key])
            pval = np.abs(np.log10(pval))
            if pvals[key] < pval:
                pvals[key] = pval

    return pvals'''

def test_igraph_partitions(parition,pheno_preval,pheno_array):
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


def test_set_partitions(partition,pheno_preval, pheno_array):
    pvals = np.zeros(len(pheno_preval))
    for subgraph in partition:
        result_dict = {}
        for vertex in subgraph:
            if pheno_array[vertex] != -1:
                if pheno_array[vertex] in result_dict:
                    result_dict[pheno_array[vertex]] += 1
                else:
                    result_dict[pheno_array[vertex]] = 1
        cls_size = len(subgraph)
        for key in result_dict:
            pval = binom_test(result_dict[key],cls_size,pheno_preval[key])
            pval = np.abs(np.log10(pval))
            if pvals[key] < pval:
                pvals[key] = pval

    return pvals

def test_purity_set(partition_list,class_list,node_count):
    total_sum = 0
    for cluster in partition_list:
        max = -1
        for class_set in class_list:
            pass

def calculate_sim(partition_list,class_list):
    sims_list = []
    for cluster in partition_list:
        temp_sim = []
        for node_set in class_list:
            temp_sim.append(len(cluster.intersection(node_set)))
        sims_list.append(temp_sim)
    return sims_list

def calculate_purity(sim_list,total_count):
    total_sum = 0
    for overlap_list in sim_list:
        total_sum += max(overlap_list)
    return total_sum/total_count

def calculate_entropy(set_list,total_count):
    return -1*np.sum([(len(cluster_set)/total_count)*np.log(len(cluster_set)/total_count) for cluster_set in set_list])


def calculate_NMI(sim_list,partition_list,class_list,total_count):
    '''mutual_info = 0
    for k in range(len(partition_list)):
        for j in range(len(class_list)):

            mutual_info += (sim_list[k][j]/total_count)*(np.log(total_count*sim_list[k][j])-np.log(len(partition_list[k])*len(class_list[j])))
    #'''
    mutual_info = np.sum([(sim_list[k][j]/total_count)*np.log(total_count*sim_list[k][j]/(len(partition_list[k])*len(class_list[j]))) for k in range(len(partition_list)) for j in range(len(class_list)) if sim_list[k][j] != 0])
    return 2*mutual_info/(calculate_entropy(partition_list,total_count)+calculate_entropy(class_list,total_count))


def test_mcl_partitions(partition,pheno_preval,pheno_array):
    pvals = np.zeros(len(pheno_preval))
    for subgraph in partition:
        result_dict = {}
        for vertext in subgraph:
            if pheno_array[vertext] != -1:
                if pheno_array[vertext] in result_dict:
                    result_dict[pheno_array[vertext]] += 1
                else:
                    result_dict[pheno_array[vertext]] = 1
        cls_size = len(subgraph)
        for key in result_dict:
            pval = binom_test(result_dict[key], cls_size, pheno_preval[key])
            pval = np.abs(np.log10(pval))
            if pvals[key] < pval:
                pvals[key] = pval
    return pvals

#def test_cluster(member_list)