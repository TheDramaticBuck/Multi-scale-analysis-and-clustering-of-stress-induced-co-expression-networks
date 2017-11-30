/* ----------------------------------------------------------------------
	
	Adapted by:

 	Dr. Nuno R. Nene, University of Cambridge, 2017

 	Email: nunonene@gmail.com


 	For specific details on the underlying methods and the datasets see:
 
        https://arxiv.org/abs/1703.02872 


	Based on/Adapted from C++ files provided with:
   
	http://www.bioconductor.org/packages/release/bioc/html/BHC.html

   
------------------------------------------------------------------------- */


#ifndef NODE_H
#define NODE_H

#include "header.h"

#include "SEGPdist.h"


class Node
{
 public:

  	Node(); 

  	static Node CreateNode(SEGPdist& segp, const int arg_nodeID);
  
  	int GetNodeID() const;

  	double GetLowerBoundLogEvidence() const;

	double GetLogDK() const;

	double GetMergePrior() const;

	


 private:

  	

  	int nodeID;
 
  	double log_d_k;
 
  	double LogEvidence;

  	double lowerBoundLogEvidence;

  	double mergePrior;
};
#endif // NODE_H
