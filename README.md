# Local IBD Clustering

This repository is created provide public access to the software material for our paper ["Selecting clustering algorithms for IBD mapping"](https://www.biorxiv.org/content/10.1101/2021.08.11.456036v1).

Supplementary Figures and Methods are available in the sca_supmat.pdf file. 

- ## run_sim.py

    In this file, we have added the code to conduct a single simulation and analysis as an example.

- ## igraph_reader.py

    The script for reading weighted edge list files into iGraph for further analysis and visualization

- ## igraph_pval.py

    Given a graph clustering from iGraph and original ground truth cluster array, this file calculate the statistical power of the clustering algorithm in a simulation.

- ## hcs2.py

    In this file, we have update the HCS algorithm implementation from this [repo](https://github.com/53RT/Highly-Connected-Subgraphs-Clustering-HCS) to support the latest version of NetworkX python library. We have also added the ability to analyze disconnected graphs to the algorithm which is essential for IBD clustering applications.

- ## clustering_tests.py

    This file includes the script that enable the calculation of clustering scores such as NMI and statistical power in the real_data_analyzer.py file.

- ## real_data_analyzer.py

    This is the main file that include most of the logic of the simulations. This file being with a set of functions that are in charge of converting the various graph formats such as `tuple_2_vertex` function that converts the output of MCL clustering library to iGraph partition objects.

  - ***Local IBD Graph Class***: This class reads a single local IBD graph file and analyzes it using the selection clustering algorithms in the `partition` function. The algorithms currently include: *MCL, Infomap, Louvain, Leiden, and HCS*. It records the time it takes each algorithm to analyze the graph. `calculate_stats` function is in charge of recording clustering metrics for each algorithm based on its performance. When no ground truth is available, this class only calculate structural metrics.

  - ***Local IBD Generator Class***: This class generetes local IBD graph files using the benchmark procedure layed out in our paper.