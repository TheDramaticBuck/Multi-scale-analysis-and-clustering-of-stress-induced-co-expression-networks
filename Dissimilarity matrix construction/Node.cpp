/* ----------------------------------------------------------------------
	
	Adapted by:

 	Dr. Nuno R. Nene, University of Cambridge, 2017

 	Email: nunonene@gmail.com


 	For specific details on the underlying methods and the datasets see:
 
        https://arxiv.org/abs/1703.02872 


	Based on/Adapted from C++ files provided with:
   
	http://www.bioconductor.org/packages/release/bioc/html/BHC.html

   
------------------------------------------------------------------------- */


#include "Node.h"
#include <limits>
extern "C" {
#include <math.h>
#include "mex.h"
}



/* ----------------------------------------------------------------------
    Default constructor.
------------------------------------------------------------------------- */

Node::Node(){}

/* ----------------------------------------------------------------------
   Factory method to instantiate this class. It creates a node consisting
   of a single data item (as specified in the argument).
------------------------------------------------------------------------- */

Node Node::CreateNode(SEGPdist& segp, const int arg_nodeID)
{
  
	Node thisNode=Node();
  
	double lengthScale, noiseFreeScale, noiseSigma;

  
	// Set all the variables to their default values

	thisNode.nodeID=arg_nodeID;

  	thisNode.log_d_k=log(dirichletProcessParameter);

  	
  	thisNode.mergePrior=1.0;
  

	thisNode.LogEvidence=-numeric_limits<double>::infinity();

  	// Find the optimised hyperparameters for this node and compute the
  	// overall log-evidence estimate

	vector<int> node=vector<int> (1);

	node[0]=arg_nodeID;

  	thisNode.lowerBoundLogEvidence
    	= segp.PairLogEvidence(node,lengthScale,noiseFreeScale,noiseSigma);

  return thisNode;
}

/* ---------------------------------------------------------------------- */

int Node::GetNodeID() const
{
  return nodeID;
}

/* ---------------------------------------------------------------------- */


double Node::GetLowerBoundLogEvidence() const
{
  return lowerBoundLogEvidence;
}

/* ---------------------------------------------------------------------- */

double Node::GetLogDK() const
{
  return log_d_k;
}
/* ---------------------------------------------------------------------- */


double Node::GetMergePrior() const
{
  return mergePrior;
}


