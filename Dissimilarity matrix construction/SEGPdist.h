/* ----------------------------------------------------------------------
	
	Adapted by:

 	Dr. Nuno R. Nene, University of Cambridge, 2017

 	Email: nunonene@gmail.com


 	For specific details on the underlying methods and the datasets see:
 
        https://arxiv.org/abs/1703.02872 


	Based on/Adapted from C++ files provided with:
   
	http://www.bioconductor.org/packages/release/bioc/html/BHC.html

   
------------------------------------------------------------------------- */


#ifndef SEGPdist_H
#define SEGPdist_H

#include <iostream>
#include <string>
#include <vector>

#include "header.h"

#include "BlockCovarianceMatrix.h"



class SEGPdist
{

	

	public:

	int noise_mode; // indicates whether there is precalculated fixed noise; for the paper 0

  	int reps;//the number of replicates per observation; for the paper =1
  
  	vector<double>  noiseData; // noise values (computed from estimators)

  	int nDataItems; // number of genes

  	int nFeatures; // number of timepoints

  	int nTimePoints;

	
  
	SEGPdist(); // constructor
  
	SEGPdist(const vector<vector<double> >& inputData);
  
	double PairLogEvidence(vector<int>& itemIndex,
				  double& lengthScale,
				  double& NoiseFree,
				  double& Noise);

	int Get_nDataItems();

	int Get_nFeatures();

	void ReadInTimePoints(vector<double> timePoints_copy);

	void ReadInNoise(vector<double> noise_copy);

	void SetNoiseMode(int mode);

	void SetReps(int num_reps);

	int GetNoiseMode();

	private:

  	
 
  	double log_d_k;
 
  	double clusterLogEvidence;

  	double lowerBoundLogEvidence;


	protected:

	vector<vector<double> > data; // nGenes * nTimePoints

  	vector<double> timePoints;

  	double dataRange; // the range of the y data ymax-ymin


	double ComputeLogDeterminant(double* choleskyMatrix, int nVariables);


	double ComputeLogEvidence(BlockCovarianceMatrix blockMatrix,
				    vector<double> data);


	void OptimiseHyperparameters(const vector<double>& yValues,
			       double& lengthScale,
			       double& noiseFreeScale,
			       double& noiseSigma);
	

  	BlockCovarianceMatrix AddNoiseToCovarianceFunction(BlockCovarianceMatrix blockMatrix, 
						     double noiseSigma);


	double ComputeGradient(const BlockCovarianceMatrix& inverseCovarianceFunction,const BlockCovarianceMatrix& covarianceDerivative,const vector<double>& alpha);


	double ComputeNoiseGradient(const BlockCovarianceMatrix& inverseCovarianceFunction,const vector<double>& alpha,double noiseSigma);


	double ComputeMaximisedLogEvidence(vector<double> yValues,
				     double& lengthScale,
				     double& noiseFreeScale,
				     double& noiseSigma);


  	BlockCovarianceMatrix SquareExponentialCovarianceFunction(double lengthScale,
							    int blockSize,
							    double noiseFreeScale);

  
	BlockCovarianceMatrix SquareExponential_lengthDerivative(double lengthScale,
							   int blockSize,
							   double noiseFreeScale);

	double ComputeLogEvidenceFromHyperparameters(const vector<double>& yValues,
					     const int blockSize,
					     const vector<double>& params);


  	void ComputeGradientsFromHyperparameters(const vector<double>& yValues,
					   const int blockSize,
					   const vector<double>& params,
					   vector<double>& grad);


  	void ComputeLogEvidenceAndGradientsFromHyperparameters(const vector<double>& yValues,
						      const int blockSize,
						      const vector<double>& params,
						      double& logEv,
						      vector<double>& grad);

  	void ImposeConstraintsOnHyperparameters(vector<double>& params);


	void LineSearch(vector<double>& xold, const double fold, vector<double>& g,
		  vector<double>& p, vector<double>& x, double& f, const double stpmax,
		  bool& check, const int blockSize, const vector<double>& yValues);


  	void DFPMaximise(vector<double>& p, const vector<int>& fix, const double gtol,
		   double& fret, const int blockSize, const vector<double>& yValues);
};
#endif // SEGPdist_H
