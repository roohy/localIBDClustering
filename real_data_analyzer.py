import numpy as np
from numpy.core.arrayprint import dtype_is_implied
import networkx as nx
import pandas as pd
import time
from hcs2 import improved_labelled_HCS
import louvain
import leidenalg
import seaborn as sns
from matplotlib import pyplot as plt
from networkx.algorithms import approximation
import markov_clustering as mcl
import igraph
from clustering_tests import test_igraph_partitions,calculate_sim,calculate_NMI,calculate_purity
from scipy.stats import nbinom,geom,poisson,expon,gamma,bernoulli,norm,binom,binom_test
import pickle as pkl
from sklearn.metrics.cluster import adjusted_mutual_info_score

THRESHOLD = 10
PVAL = 0.0000000001

#changing MCL outout to igraph clustering format
def tuple_2_vertex(mcl_clusters,graph):
    cluster_count = 0 
    membership_array = np.zeros((len(graph.vs)),dtype=np.int)
    for i in range(len(mcl_clusters)):
        for vertex in mcl_clusters[i]:
            membership_array[vertex] = cluster_count
        cluster_count += 1
    return igraph.clustering.VertexClustering(graph,membership_array)

#^you can do this since sparse array generation doesn't mess with the order of the nodes, making it possible to recover the names based on order alone
#changing a set partition to a membership array format. Membership arrays can be easily converted into igraph clustering
def set_2_membership_array(partition,nodeCount):
    membershipArray = np.zeros((nodeCount))
    clusterCounter = 0 
    for cluster in partition:
        for node in cluster:
            membershipArray[node] = clusterCounter
        clusterCounter += 1
    return membershipArray

#igraph vertex clustering to set partition 
def vertex_to_set(vertex_cluster):
    parition = []
    for subgraph in vertex_cluster.subgraphs():
        subgraph_set = set()
        for vertex in subgraph.vs:
            subgraph_set.add(int(vertex['name']))
        if len(subgraph_set) > 0:
            parition.append(subgraph_set)
    return parition

#igraph vertex clustering to membership array
def vertex_to_list(vertex_cluster,node_count):
    membership_array = np.zeros((node_count),dtype=int)
    counter = 1
    for subgraph in vertex_cluster.subgraphs():
        for vertex in subgraph.vs:
            membership_array[int(vertex['name'])] = counter
        counter += 1 
    return membership_array
#deprecated    
def update_membership_array(membership_array, igraph_object):
    updated_membership_array = []
    for vertex in igraph_object.vs:
        updated_membership_array.append(membership_array[int(vertex['name'])])
    return np.array(updated_membership_array)

#Main class that reads a graph file and runs the clustering algorithms over it.
class LocalIBDGraph(object):
    clustering_methods = ['mcl1.5','mcl2','mcl3','mcl5','louvain','leiden','infomap','fake_infomap']#,'HCS'
    
    #pass a True value for simulated if a ground truth file is available. The address to true info should be passed in pheno_path
    def __init__(self,addr,simulated=False,pheno_path=None):
        super().__init__()
        self.igg = igraph.Graph.Read_Ncol(addr,directed=False)
        self.nxg = nx.read_edgelist(addr)
        self.nodeCount = len(self.nxg)
        if simulated:
            self.truth = pkl.load(open(pheno_path,'rb'))
        else:
            self.truth = None


    def make_real_clustering(self):
        # self.membership_array = update_membership_array(self.truth['membership_array'],self.igg)
        updated_membership_array = []
        for vertex in self.igg.vs:
            updated_membership_array.append(self.truth['membership_array'][int(vertex['name'])])
        self.membership_array = np.array(updated_membership_array)
        self.real_membership_array = self.truth['membership_array']
    
    def partition(self):
        self.clustering_results = []
        
        
        sparseMat =  nx.to_scipy_sparse_matrix(self.nxg)
        self.mclPartitionRes = []
        self.clustering_stats = {}
        for i in [1.5,2,3,5]:
            tic = time.perf_counter()
            self.clustering_results.append(tuple_2_vertex(mcl.get_clusters(mcl.run_mcl(sparseMat,inflation=i)),self.igg))
            toc = time.perf_counter()
            self.clustering_stats['mcl{}_time'.format(i)] = toc-tic
        tic = time.perf_counter()
        self.clustering_results.append(louvain.find_partition(self.igg,louvain.ModularityVertexPartition))
        toc = time.perf_counter()
        self.clustering_stats['louvain_time'] = toc-tic
        self.clustering_results.append(leidenalg.find_partition(self.igg,leidenalg.ModularityVertexPartition))
        tic = time.perf_counter()
        self.clustering_stats['leiden_time'] = tic-toc
        self.clustering_results.append(self.igg.community_infomap())
        toc = time.perf_counter()
        self.clustering_stats['infomap_time'] = toc-tic
        # if self.truth is None:
        #     self.clustering_results.append(self.igg.components())
        
        #hcs_mem_array = improved_labelled_HCS(self.nxg)
        #hcs_mem_array = update_membership_array(hcs_mem_array,self.igg)
        #self.clustering_results.append(igraph.clustering.VertexClustering(self.igg,hcs_mem_array))
        tic = time.perf_counter()
        #self.clustering_stats['HCS_time'] = tic-toc
        self.clustering_results.append(igraph.clustering.VertexClustering(self.igg,self.randomize(self.clustering_results[-1])))
        #self.clustering_results.append(self.igg.components()) #Decide later 
        self.clustering_stats['node_count'] = self.nodeCount
        self.clustering_stats['edge_count'] = self.nxg.number_of_edges()
        
        if self.truth is not None:
            self.make_real_clustering()

            # self.clustering_results.append(igraph.clustering.VertexClustering(self.igg,
            #      improved_labelled_HCS(self.nxg,self.membership_array.shape[0])))
            self.real_clustering = igraph.clustering.VertexClustering(self.igg, self.membership_array.astype(int))
            self.clustering_stats['real_modularity'] = self.real_clustering.modularity
        
        for i in range(len(self.clustering_methods)):
            self.clustering_stats[self.clustering_methods[i]+'_modularity']=self.clustering_results[i].modularity
    
    
    def calculate_stats(self):
        if self.truth is not None:
            self.get_pvals()
            self.calculate_gt_stats()
        for i in range(len(self.clustering_methods)): 
            self.clustering_stats.update(self.calculate_stats_igraph(self.clustering_results[i],self.clustering_methods[i]))
        return self.clustering_stats
    
    
    
    
    
    def calculate_stats_igraph(self,partition,clusteringMethodName):
        resultDict = {}
        nodeCount = []

        totalConnectivity = 0 
        highlyConnectedClusters = 0
        largeClusterCount = 0 
        clusterCount = 0 
        outsideEdgeCount = 0
        missingEdgeCount = 0
        expectedEdgeCount = 0
        trueEdgeCount = 0 
        for subgraph in partition.subgraphs():
            clusterCount += 1 

            clusterNodeCount = len(subgraph.vs)
            
            vertexSet = set()
            edgeCount = 0
            if clusterNodeCount == 1:
                outsideEdgeCount += self.nxg.degree(subgraph.vs[0]['name'])
                continue
            nodeCount.append(clusterNodeCount)
            for source in subgraph.vs:
                vertexSet.add(source['name'])
            for source in vertexSet:
                nodeConnectivity = 0
                nodeDegree = self.nxg.degree(source)
                for target in vertexSet:
                    if self.nxg.has_edge(source,target):
                        nodeConnectivity+= 1
                        edgeCount += 1
                        nodeDegree -= 1
                
                if nodeConnectivity >= clusterNodeCount/2:
                    totalConnectivity += 1
                outsideEdgeCount += nodeDegree
            if clusterNodeCount > 5:
                largeClusterCount += 1
                if clusterNodeCount*(clusterNodeCount-1)/2 < edgeCount:
                    highlyConnectedClusters += 1
            expectedEdgeCount += clusterNodeCount*(clusterNodeCount-1)    
            missingEdgeCount += (clusterNodeCount*(clusterNodeCount-1)) - edgeCount
            trueEdgeCount += edgeCount
        nodeCount = np.array(nodeCount)

        degreeCount = nodeCount*(nodeCount-1)
        
        resultDict[clusteringMethodName+'_cls_count'] = clusterCount
        resultDict[clusteringMethodName+'_node_count'] = nodeCount.sum()
        resultDict[clusteringMethodName+'_degree_count'] = degreeCount.sum()
        resultDict[clusteringMethodName+'_node_per_cluster'] = nodeCount.mean()
        resultDict[clusteringMethodName+'_node_per_cluster_std'] = nodeCount.std()
        resultDict[clusteringMethodName+'_total_connectivity'] = totalConnectivity
        resultDict[clusteringMethodName+'_connectivity'] = totalConnectivity/nodeCount.sum()
        resultDict[clusteringMethodName+'_outside_edge_count'] = outsideEdgeCount/2
        resultDict[clusteringMethodName+'_outside_edge_ratio'] = outsideEdgeCount/(2*self.nxg.number_of_edges())
        resultDict[clusteringMethodName+'_highly_connected_count'] = highlyConnectedClusters
        resultDict[clusteringMethodName+'_highly_connected_ratio'] = 0 if largeClusterCount == 0 else highlyConnectedClusters/largeClusterCount
        resultDict[clusteringMethodName+'_large_cluster_count'] = largeClusterCount
        resultDict[clusteringMethodName+'_missing_edge_count'] = missingEdgeCount/2
        resultDict[clusteringMethodName+'_missing_edge_ratio'] = missingEdgeCount/expectedEdgeCount
        resultDict[clusteringMethodName+'_missing_edge_over_available_edges'] = missingEdgeCount/degreeCount.sum()
        resultDict[clusteringMethodName+'_coverage'] = trueEdgeCount/(2*self.nxg.number_of_edges())
        # resultDict[clusteringMethodName+'_']
        return resultDict
    
    
    def get_pvals(self):
        clustering_pvals = {}
        for i in range(len(LocalIBDGraph.clustering_methods)):
            clustering_pvals[LocalIBDGraph.clustering_methods[i]] = test_igraph_partitions(self.clustering_results[i],
                self.truth['pheno_preval'],self.truth['phenotyep_array'])
        max_pval = np.sum(test_igraph_partitions(self.real_clustering,self.truth['pheno_preval'],self.truth['phenotyep_array']))
        for i in range(len(LocalIBDGraph.clustering_methods)):
            pvals = clustering_pvals[LocalIBDGraph.clustering_methods[i]]
            if max_pval == 0:
                self.clustering_stats[LocalIBDGraph.clustering_methods[i]+'_sig_score'] = -1
                self.clustering_stats[LocalIBDGraph.clustering_methods[i]+'_bc_sig_score'] = -1
            else:
                self.clustering_stats[LocalIBDGraph.clustering_methods[i]+'_sig_score'] = np.sum(pvals)/max_pval
                self.clustering_stats[LocalIBDGraph.clustering_methods[i]+'_bc_sig_score'] = np.sum(pvals[pvals>=5])/max_pval
    
    
    
    def calculate_gt_stats(self):
        clustering_set_results = []
        for cls_result in self.clustering_results:
            clustering_set_results.append(vertex_to_set(cls_result))
        # self.cluster_sim_results = [] 
        for i in range(len(clustering_set_results)):
            sim_result = calculate_sim(clustering_set_results[i],self.truth['cls_dict'])
            self.clustering_stats[self.clustering_methods[i]+'_purity'] = calculate_purity(sim_result,self.truth['total_count'])
            self.clustering_stats[self.clustering_methods[i]+'_NMI'] = calculate_NMI(sim_result,clustering_set_results[i],
                self.truth['cls_dict'],self.truth['total_count'])
            temp_membership_array = vertex_to_list(self.clustering_results[i],self.real_membership_array.shape[0])
            self.clustering_stats[self.clustering_methods[i]+'_AMI'] = adjusted_mutual_info_score(self.real_membership_array,
                temp_membership_array)
    
    
    def randomize(self,partition):
        clusterCount = 0
        nodeCounts = []
        for subgraph in partition.subgraphs():
            clusterCount += 1 

            clusterNodeCount = len(subgraph.vs)
            nodeCounts.append(clusterNodeCount)
        membership_array = np.zeros((self.nodeCount),dtype=int)
        permutation = np.random.permutation(self.nodeCount)
        head = 0 
        clusterCounter = 1
        for subgraphNodeCount in nodeCounts:
            membership_array[permutation[head:head+subgraphNodeCount]] = clusterCounter
            clusterCounter += 1
            head += subgraphNodeCount
        return membership_array

class LocalIBDGenerator(object):
    def __init__(self,cluster_distro):
        self.cluster_distro = cluster_distro    
        self.create_graphs()
    def create_graphs(self):
        inds = np.sum(self.cluster_distro)
        cls_count = self.cluster_distro.shape[0]
        self.membership_array = np.zeros((inds))
        index = 0
        sample_counter = 0
        self.cls_dict = []
        for graph_size in self.cluster_distro:
            temp_array = set()
            for i in range(graph_size):
                temp_array.add(sample_counter)
                self.membership_array[sample_counter] = index
                sample_counter += 1
            self.cls_dict.append(temp_array)
            index += 1
        self.total_count = inds
        return self.membership_array, self.cls_dict
    def generate_file(self, fp_rate, fn_rate,path):
        with open(path+'.abc','w') as output_file:
            people_count = self.membership_array.shape[0]
        
            counter = 0
            for cls_set in self.cls_dict:
                cls_list = list(cls_set)
                for i in range(len(cls_list)):
                    for j in range(i + 1, len(cls_list)):
                        if not bernoulli.rvs(fn_rate):
                            output_file.write('\t'.join([str(item) for item in [cls_list[i],cls_list[j]]])+'\n')
                        
                            counter += 1
            counter = int(fp_rate * counter)
            while counter > 0:
                chosen_ones = np.random.choice(people_count, 2, replace=False)
                if self.membership_array[chosen_ones[0]] != self.membership_array[chosen_ones[1]]:
                    output_file.write('\t'.join([str(item) for item in [chosen_ones[0],chosen_ones[1]]])+'\n')
                    counter -= 1
    def create_phenotype(self,general_mean,general_std):
        phenotype_array = np.zeros(self.membership_array.shape[0],dtype=object)
        for i in range(self.membership_array.shape[0]):
            phenotype_array[i] = []
        prob_array = np.zeros(self.membership_array.shape[0],dtype=float)
        self.pheno_preval = []
        pheno_count = 0 
        for i in range(len(self.cls_dict)):
            cls_listed = list(self.cls_dict[i])
            if len(self.cls_dict[i]) >= THRESHOLD:
                #temp_pheno_array = np.zeros(self.membership_array.shape[0],dtype=bool)
                #temp_pheno_array[:] = False
                preval = np.abs(norm.rvs(loc=general_mean,scale=general_std))
                general_count = self.membership_array.shape[0]-len(self.cls_dict[i])
                general_affected_count = binom.rvs(general_count,preval)
                prob_array[:] = 1
                prob_array[cls_listed] = 0
                prob_array = prob_array/general_count
                general_index_list = np.random.choice(self.total_count,general_affected_count,replace=False,p=prob_array)
                for item in general_index_list:
                    phenotype_array[item].append(pheno_count)
                
                affected_count = int(binom.isf(PVAL,len(self.cls_dict[i]),preval))
                index_list = np.random.choice(len(self.cls_dict[i]),affected_count,replace=False)
                for index in index_list:
                    phenotype_array[cls_listed[index]].append(pheno_count)
                    
                pheno_count += 1
                self.pheno_preval.append(preval)
                
        
        
        
        
        self.phenotype_array = phenotype_array
        
        return self.phenotype_array,self.pheno_preval
            

    def write_to_file(self,path):
        result = {'total_count':self.total_count,'cls_dict':self.cls_dict,
                    'membership_array':self.membership_array,'pheno_preval':self.pheno_preval,
                    'phenotyep_array':self.phenotype_array}
        pkl.dump(result,open(path+'pkl','wb'))
            
                
