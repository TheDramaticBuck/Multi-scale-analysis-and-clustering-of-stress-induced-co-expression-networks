/* ----------------------------------------------------------------------
	
	Adapted by:

 	Dr. Nuno R. Nene, University of Cambridge, 2017

 	Email: nunonene@gmail.com


 	For specific details on the underlying methods and the datasets see:
 
        https://arxiv.org/abs/1703.02872 


	Based on/Adapted from C++ files provided with:
   
	http://www.bioconductor.org/packages/release/bioc/html/BHC.html

   
------------------------------------------------------------------------- */


#ifndef GPdist_H
#define GPdist_H

#include "header.h"

#include "Node.h"
#include "SEGPdist.h"



class GPdist
{
 public:

  	GPdist();

  	vector<double> ComputeGPdist(SEGPdist& segp);

	double ComputeLikPair(SEGPdist& segp,
			    const Node& node1,
			    const Node& node2);
  
  
 private:


	
 
  	void ComputeGPdistPairWise(vector<double>& LogEvidence,vector<Node>& network,
		 SEGPdist& segp);

};

#endif // GPdist_H
