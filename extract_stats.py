import numpy as np
from numpy.core.arrayprint import dtype_is_implied
import networkx as nx
import igraph,sys,pickle


class LocalIBDGraph(object):
    #pass a True value for simulated if a ground truth file is available. The address to true info should be passed in pheno_path
    def __init__(self,cls_addr,edge_list):
        self.cls_list = []
        with open(cls_addr,'r') as cls_file:
            for line in cls_file:
                data = line.strip().split()
                self.cls_list.append(data)
        self.igg = igraph.Graph.Read_Ncol(edge_list,directed=False)
        self.nxg = nx.read_weighted_edgelist(edge_list)
        self.nodeCount = len(self.nxg)
        
    def make_partition(self):
        self.map = {}
        for index,vertex in enumerate(self.igg.vs):
            self.map[vertex['name']] = index
        self.membership_array = np.zeros((len(self.igg.vs)),dtype=np.int)
        for index,cluster in enumerate(self.cls_list):
            for vertex in cluster:
                self.membership_array[self.map[vertex]] = index
        self.partition = igraph.clustering.VertexClustering(self.igg,self.membership_array)
        return self.partition

    def calculate_stats_igraph(self,partition):
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
        
        resultDict['cls_count'] = clusterCount
        resultDict['node_count'] = nodeCount.sum()
        resultDict['degree_count'] = degreeCount.sum()
        resultDict['node_per_cluster'] = nodeCount.mean()
        resultDict['node_per_cluster_std'] = nodeCount.std()
        resultDict['total_connectivity'] = totalConnectivity
        resultDict['connectivity'] = totalConnectivity/nodeCount.sum()
        resultDict['outside_edge_count'] = outsideEdgeCount/2
        resultDict['outside_edge_ratio'] = outsideEdgeCount/(2*self.nxg.number_of_edges())
        resultDict['highly_connected_count'] = highlyConnectedClusters
        resultDict['highly_connected_ratio'] = 0 if largeClusterCount == 0 else highlyConnectedClusters/largeClusterCount
        resultDict['large_cluster_count'] = largeClusterCount
        resultDict['missing_edge_count'] = missingEdgeCount/2
        resultDict['missing_edge_ratio'] = missingEdgeCount/expectedEdgeCount
        resultDict['missing_edge_over_available_edges'] = missingEdgeCount/degreeCount.sum()
        resultDict['coverage'] = trueEdgeCount/(2*self.nxg.number_of_edges())
        # resultDict[clusteringMethodName+'_']
        return resultDict
    
    
if __name__ == '__main__':

    edge_list_addr = sys.argv[1]
    cls_file_addr = sys.argv[2]
    output_addr = sys.argv[3]

    libdg = LocalIBDGraph(cls_file_addr,edge_list_addr)
    partition = libdg.make_partition()
    stats = libdg.calculate_stats_igraph(partition)
    pickle.dump(stats,open(output_addr,'wb'))