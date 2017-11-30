/* ----------------------------------------------------------------------
	
	Adapted by:

 	Dr. Nuno R. Nene, University of Cambridge, 2017

 	Email: nunonene@gmail.com


 	For specific details on the underlying methods and the datasets see:
 
        https://arxiv.org/abs/1703.02872 


	Based on/Adapted from C++ files provided with:
   
	http://www.bioconductor.org/packages/release/bioc/html/BHC.html

   
------------------------------------------------------------------------- */


#include <limits>



#include "GPdist.h"

extern "C" {
#include <math.h>
#include "mex.h"
}




/* ----------------------------------------------------------------------
   Default constructor.
---------------------------------------------------------------------- */
   
GPdist::GPdist() {}


/* ----------------------------------------------------------------------
   Compute joint likelihood
------------------------------------------------------------------------- */

double GPdist::ComputeLikPair(SEGPdist& segp,
			    const Node& node1,
			    const Node& node2)
{

  	double tr1,tr2,a,b,ckt,pk,gell,num1,num2;

	vector<int> NodePairIDs=vector<int>(2);
	
	NodePairIDs[1]=node1.GetNodeID(); 

	NodePairIDs[2]=node2.GetNodeID();

  	double lengthScale, noiseFreeScale, noiseSigma;

  	// Compute the pair of nodes log evidence and lower bound log evidence
  
	tr1=log(dirichletProcessParameter)+gammaln(node1.GetMergePrior()+node2.GetMergePrior());
  
	tr2=node1.GetLogDK()+node2.GetLogDK();
 	
	a=max(tr1,tr2);
  
	b=min(tr1,tr2);
  
	ckt=a+log(1.0+exp(b-a));
  
	pk=tr1-ckt;
  
	gell=segp.PairLogEvidence(NodePairIDs,
					lengthScale,
					noiseFreeScale,
					noiseSigma);

	num1=pk+gell;
  
	num2=tr2-ckt+node1.GetLowerBoundLogEvidence()+node2.GetLowerBoundLogEvidence();

  return num1-num2;
}


/* ----------------------------------------------------------------------
   Compute GP distance matrix
---------------------------------------------------------------------- */

vector<double> GPdist::
ComputeGPdist(SEGPdist& segp)
{
  int n=segp.nDataItems;
  int i;

  vector<Node> network=vector<Node>(n);

  vector<double> LogEvidence=vector<double>(n*(n-1)/2);

  // Calculate the independent likelihood for each node on the network

#pragma omp parallel for default(shared) private(i) schedule(dynamic,1)

  for(i=0; i<n; i++) network[i]=Node::CreateNode(segp, i);

  // Compute joint likelihood

  ComputeGPdistPairWise(LogEvidence,network, segp);
  
  return LogEvidence;
}

/* ----------------------------------------------------------------------
   Calculate log of ratio of probabilities
---------------------------------------------------------------------- */

void GPdist::
ComputeGPdistPairWise(vector<double>& LogEvidence,vector<Node>& network,
		 SEGPdist& segp)
{
  int n=network.size(); // full network

  int i,j;

  // Compute all possible joint likelihood

#pragma omp parallel for default(shared) private(i,j) schedule(dynamic,1)

  for(i=1; i<n; i++)
    {
      int offset=(i-1)*i/2;

      for(j=0; j<i; j++)
	{
	  // Pair nodes i and j

	  LogEvidence[offset+j]=ComputeLikPair(segp,network[i],network[j]);

	  LogEvidence[offset+j]=LogEvidence[offset+j]-((network[i].GetLowerBoundLogEvidence())+(network[j].GetLowerBoundLogEvidence())); // log of ration between joint and separate likelihoods
	}
    }

     
}

