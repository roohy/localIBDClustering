import sys
import igraph

def vertex_to_set(vertex_cluster): #Converts iGraph partitionings to a list of clusters. Each cluster represented by a set
    parition = []
    for subgraph in vertex_cluster.subgraphs():
        subgraph_set = set()
        for vertex in subgraph.vs:
            subgraph_set.add(int(vertex['name']))
        if len(subgraph_set) > 0:
            parition.append(subgraph_set)
    return parition

if __name__ == '__main__':
    graph_addr = sys.argv[1] #Address to the graph edgelist file. Each line has two label nodes
    igg = igraph.Graph.Read_Ncol(graph_addr,directed=False) #edgelist reader function for iGraph. it generates an iGraph instance and loads the file into it.
    partition = igg.community_infomap()
    set_partition = vertex_to_set(partition)
    #...