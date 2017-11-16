
Homegeneity, Separation and Gene Ontology Overlap
=================================================

This function computes the value of linear performance functions between gene expression profiles grouped in clusters

Clustering solutions were computed through stability analysis

Clustering performance functions: Hav (average homogeneity), Sav (average separation) and 
GOTOav (average gene ontology term overlap) at each resolution, here quantified by the Markov Time

For specific details on the underlying methods, including stability analysis, relevant references and the datasets see:
 
        https://arxiv.org/abs/1703.02872 

function [Hav_MT,Sav_MT,GOTOavbp_MT,GOTOavmf_MT,GOTOavcc_MT]=HavSavGOTOav(T,probenames,datamatrix,SGDmap_bp,SGDmap_mf,SGDmap_cc,C)

Inputs:

T .......................................... Markov Times for each partition solution determined through stability analysis
probenames ................................. Identifiers for each row of datamatrix
datamatrix ................................. Row-wise normalized expression matrix: n probenames by p time-points
SGDmap_bp, SGDmap_mf and SGDmap_cc ......... Map container of Gene Ontology terms for S.cerevisiae (SGD); import with load('SC_CompleteGOterms.mat')
C .......................................... n probes by length(T) matrix of partition solutions identified by stability analysis; each column is a vector where entries are identified by their respective cluster or community at the respective Markov time

Outputs:

Hav_MT ..................................... Average homegeneity for each Markov time in T
Sav_MT ..................................... Average separation for each Markov time in T
GOTObp_MT .................................. Average GOTO (bp) for each Markov time in T
GOTOmf_MT .................................. Average GOTO (mf) for each Markov time in T
GOTOcc_MT .................................. Average GOTO (cc) for each Markov time in T
