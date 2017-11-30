# Multi-scale-analysis-and-clustering-of-stress-induced-co-expression-networks
==============================================================================

Code for several stages of a protocol for clustering analysis and integration of multiple dataset multi-scale organisation under a graph-theoretical paradigm; described in https://arxiv.org/abs/1703.02872.

PaperWithSupplFigs folder
-------------------------
Paper with all supplementary figures included. Extended version of arxiv: https://arxiv.org/abs/1703.02872. 


Dissimilarity matrix and co-expression network construction folder
------------------------------------------------------------------
Contains code for the construction of dissimilarity matrices according to the functions outlined in the paper and an implementation of the RMST algorithm in C++ (see paper for details).


Clustering performance evaluation folder
----------------------------------------
Routines in C++ for evaluation of clustering performance (LOOCVav)(available soon).
Routines in Matlab for linear performance functions: Homogeneity, Separation and Gene Ontology Overlap.

AMI folder
----------

Matlab files for computing the distance in partition space (1-AMI) between clustering solutions determined by stability analysis.

Run_CLICK_R_Linux
-----------------

Further results generated through the CLICK algorithm (not provided in the paper): http://acgt.cs.tau.ac.il/expander/.
Tested only in Linux.
