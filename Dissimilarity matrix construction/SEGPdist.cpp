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

#include "SEGPdist.h"

#include "BlockCovarianceMatrix.h"

extern "C" {
#include <math.h>
#include "mex.h"
}


/* ---------------------------------------------------------------------- */
   
SEGPdist::SEGPdist() {}


/* ---------------------------------------------------------------------- */


SEGPdist::SEGPdist(const vector<vector<double> >& inputData)
{
  
  data = inputData;

  nDataItems  = data.size();

  nFeatures   = data[0].size();

  nTimePoints = nFeatures;

  
}


/* ---------------------------------------------------------------------- */

int SEGPdist::Get_nDataItems()
{
  return nDataItems; // number of genes
}

/* ---------------------------------------------------------------------- */

int SEGPdist::Get_nFeatures()
{
  return nFeatures; // number of timepoints
}



/* ---------------------------------------------------------------------- */

void SEGPdist::ReadInTimePoints(vector<double> timePoints_copy)
{

  int i;
  
  for (i=0; i<nTimePoints; i++)
  {
    timePoints.push_back(timePoints_copy[i]);
  }
}

/* ---------------------------------------------------------------------- */

void SEGPdist::ReadInNoise(vector<double> noise_copy)
{
  noiseData = noise_copy;
}

/* ---------------------------------------------------------------------- */

void SEGPdist::SetNoiseMode(int mode)
{
  noise_mode = mode;

}

/* ---------------------------------------------------------------------- */

void SEGPdist::SetReps(int num_reps)
{
  reps = num_reps;
 
}


/* ---------------------------------------------------------------------- */

int SEGPdist::GetNoiseMode()
{
  return(noise_mode);
}

/* ---------------------------------------------------------------------- */

double SEGPdist::ComputeLogEvidence(BlockCovarianceMatrix blockMatrix,
					     vector<double> data)
{
  //DECLARATIONS

  double logEvidence;

  const double PI=3.14159265358979324;

  logEvidence  = -0.5 * blockMatrix.ComputeMatrixDeterminant();

  logEvidence -= 0.5 * nTimePoints * blockMatrix.blockSize * log(2*PI);

  //invert the BlockMatrix
  blockMatrix.InvertMatrix();
  //compute the likelihood term
  logEvidence -= 0.5*blockMatrix.ComputeLogLikelihoodProduct(data);

  //IT MIGHT BE SENSIBLE TO FORBID logEv=inf HERE (RESULT OF A SINGULAR MATRIX)

  //ARE THERE ANY DANGERS TO DOING THIS??
  if (logEvidence==numeric_limits<double>::infinity())
    logEvidence = -numeric_limits<double>::infinity();//-ve inf gives us ev=0, which will be rejected

  return(logEvidence);
}

/* ---------------------------------------------------------------------- */

double SEGPdist::
ComputeGradient(const BlockCovarianceMatrix& inverseCovarianceFunction,
		const BlockCovarianceMatrix& covarianceDerivative,
		const vector<double>& alpha)
{
  // At first glance, it would appear to be quicker to explicitly build the
  // rows as you go along; but doing it this way cannot exploit the STL
  // and so it turns out quicker to just use GetRow() as below - R.D.

  // Declarations
  size_t k;
  int j;
  double gradient=0.0;
  const int alpha_size=alpha.size();
  vector<double> vec1=vector<double>(alpha_size);
  vector<double> vec2=vector<double>(alpha_size);
  double alpha_dot_vec2, vec1_dot_vec2;
  int block_j1, block_j2, block_j1_counter, block_j2_counter;
  vector<double>::const_iterator inIt1, inIt2, ub;
  vector<double>::iterator vecIt1, vecIt2, outItEnd1, outItEnd2;
  
  // Compute the gradient
  block_j1_counter = block_j2_counter = block_j1 = block_j2 = 0;
  for (j=0; j<alpha_size; j++)
  {
    //vec1 = inverseCovarianceFunction.GetRow(j);
    //vec2 = covarianceDerivative.GetRow(j);
    // ^ we merge these two operations into one (convoluted but quick):
    ///////////////////////////////////////////////////////////////////////
    inIt1 = inverseCovarianceFunction.noiseFreeCoeff[block_j1].begin();
    vecIt1 = vec1.begin();
    inIt2 = covarianceDerivative.noiseFreeCoeff[block_j2].begin();
    vecIt2 = vec2.begin();
    ub = inverseCovarianceFunction.noiseFreeCoeff[block_j1].end();
    while (inIt1 != ub)
      {
	outItEnd1 = vecIt1 + inverseCovarianceFunction.blockSize;
	outItEnd2 = vecIt2 + covarianceDerivative.blockSize;
	fill(vecIt1, outItEnd1, *inIt1++);
	fill(vecIt2, outItEnd2, *inIt2++);
	vecIt1 = outItEnd1;
	vecIt2 = outItEnd2;
      }
    vec1[j] *= 1.0 + inverseCovarianceFunction.noisyCoeff[block_j1];
    vec2[j] *= 1.0 + covarianceDerivative.noisyCoeff[block_j2];
    // block_j1 = j / inverseCovarianceFunction.blockSize
    // block_j2 = j / covarianceDerivative.blockSize
    // ^ we can do this without doing the expensive divisions:
    if(++block_j1_counter==inverseCovarianceFunction.blockSize)
      {
	block_j1_counter=0; ++block_j1;
      }
    if(++block_j2_counter==covarianceDerivative.blockSize)
      {
	block_j2_counter=0; ++block_j2;
      }
    /////////////////////////////////////////////////////////////////////
    
    // Compute alpha.dK/dtheta and K^-1.dK/dtheta
    alpha_dot_vec2 = vec1_dot_vec2 = 0;
    for (k=0; k<alpha.size(); k++)
    {
      alpha_dot_vec2 += alpha[k] * vec2[k];
      vec1_dot_vec2 += vec1[k] * vec2[k];
    }
    
    // This element's contribution to the trace
    gradient += alpha[j] * alpha_dot_vec2 - vec1_dot_vec2;
  }

  return gradient * 0.5;
}

/* ---------------------------------------------------------------------- */

double SEGPdist::
ComputeNoiseGradient(const BlockCovarianceMatrix& inverseCovarianceFunction,
		     const vector<double>& alpha,
		     double noiseSigma)
{
  //DECLARATIONS
  unsigned int i;
  double gradient=0, currentElement;
  
  //COMPUTE THE GRADIENT
  for (i=0; i<alpha.size(); i++)
  {
    currentElement = alpha[i] * alpha[i];
    currentElement -= inverseCovarianceFunction.GetElement(i,i);
    gradient += currentElement;
  }
  //was: gradient * noiseSigma;
  return gradient * 0.5;
}


/* ---------------------------------------------------------------------- */

BlockCovarianceMatrix SEGPdist::
AddNoiseToCovarianceFunction(BlockCovarianceMatrix blockMatrix, double noiseSigma)
{
  // Declarations
  int i;
  double sigmaSquared;
  
  // Add sigma^2 (variange) to the diagonal elements
  sigmaSquared           = pow(noiseSigma, 2);
  for (i=0; i<blockMatrix.nRank; i++)
    // normalisation because of the way we construct a block matrix
    blockMatrix.noisyCoeff[i] = sigmaSquared / blockMatrix.noiseFreeCoeff[i][i];

  return(blockMatrix);
}



/* ---------------------------------------------------------------------- */

void SEGPdist::
LineSearch(vector<double>& xold, const double fold, vector<double>& g,
	   vector<double>& p, vector<double>& x, double& f, const double stpmax,
	   bool& check, const int blockSize, const vector<double>& yValues)
{
  const double ALF=1.0e-3, TOLX=numeric_limits<double>::epsilon();
  int i;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
  double rhs1,rhs2,slope,sum,temp,test,tmplam;

  int n=xold.size();
  check=false;
  sum=0.0;
  for(i=0; i<n; i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if(sum > stpmax)
    {
      temp=stpmax/sum;
      for(i=0; i<n; i++) p[i] *= temp;
    }
  slope=0.0;
  for(i=0; i<n; i++) slope += g[i]*p[i];
  //if(slope >= 0.0) cout << "Roundoff problem in line_search: " << slope << endl;
  test=0.0;
  for(i=0; i<n; i++)
    {
      temp=fabs(p[i])/MAX(fabs(xold[i]),1.0);
      if(temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for(;;)
    {
      for(i=0; i<n; i++) x[i]=xold[i] + alam*p[i];
      ImposeConstraintsOnHyperparameters(x);
      f=ComputeLogEvidenceFromHyperparameters(yValues,blockSize,x);
      f=-f;
      if(alam < alamin)
	{
	  for(i=0; i<n; i++) x[i]=xold[i];
	  check=true;
	  return;
	}else if(f <= fold + ALF*alam*slope) return;
      else
	{
	  if(alam == 1.0)
	    tmplam = -slope/(2.0*(f-fold*slope));
	  else
	    {
	      rhs1=f-fold-alam*slope;
	      rhs2=f2-fold-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if(a == 0.0) tmplam = -slope/(2.0*b);
	      else
		{
		  disc=b*b-3.0*a*slope;
		  if(disc < 0.0) tmplam=0.5*alam;
		  else if(b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
		  else tmplam=-slope/(b+sqrt(disc));
		}
	      if(tmplam>0.5*alam)
		tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2=f;
      alam=MAX(tmplam,0.1*alam);
    }
}

/* ----------------------------------------------------------------------
   Uses the Broyden-Fletcher-Goldfarb-Shanno (BFGS) variant of the
   Davidon-Fletcher-Powell (DFP) optimisation method (quasi-Newton) to
   perform function maximisation.

   This function can only be called from a class that inherits this class,
   and which must define the functions:
      - ComputeLogEvidenceAndGradientsFromHyperparameters()
      - ComputeGradientsFromHyperparameters()
      - ComputeLogEvidenceFromHyperparameters()
      - ImposeConstraintsOnHyperparameters()

   Adapted from Numerical Recipes in C++.
------------------------------------------------------------------------- */

void SEGPdist::
DFPMaximise(vector<double>& p, const vector<int>& fix, const double gtol, double& fret, const int blockSize, const vector<double>& yValues)
{
  const int ITMAX=100;
  const double EPS=numeric_limits<double>::epsilon();
  const double TOLX=4*EPS, STPMX=100.0;
  bool check, fixedparam=false;
  int i,its,j;
  double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

  int n=p.size();
  vector<double> dg(n), g(n), hdg(n), pnew(n), xi(n);
  vector<vector<double> > hessin=vector<vector<double> >(n,vector<double>(n,0.0));

  for(i=0; i<n; i++) if(fix[i]) {fixedparam=true; break;}
  
  ComputeLogEvidenceAndGradientsFromHyperparameters(yValues,blockSize,
						    p,fp,g);
  fp=-fp;
  for(i=0; i<n; i++) g[i]=-g[i];
  if(fixedparam) for(i=0; i<n; i++) if(fix[i]) g[i]=0;
  for(i=0; i<n; i++)
    {
      hessin[i][i]=1.0; // identity
      xi[i] = -g[i];
      sum += p[i]*p[i];
    }
  stpmax=STPMX*MAX(sqrt(sum),double(n));
  for(its=0; its<ITMAX; its++)
    {
      LineSearch(p,fp,g,xi,pnew,fret,stpmax,check,blockSize,yValues);

      fp=fret;
      for(i=0; i<n; i++)
	{
	  xi[i]=pnew[i]-p[i];
	  p[i]=pnew[i];
	}
      test=0.0;
      for(i=0; i<n; i++)
	{
	  temp=fabs(xi[i])/MAX(fabs(p[i]),1.0);
	  if(temp>test)test=temp;
	}
      if(test < TOLX) return;
      for(i=0; i<n; i++) dg[i]=g[i];
      ComputeGradientsFromHyperparameters(yValues,blockSize,
					  p,g);
      for(i=0; i<n; i++) g[i]=-g[i];
      if(fixedparam) for(i=0; i<n; i++) if(fix[i]) g[i]=0;
      test=0.0;
      den=MAX(fret,1.0);
      for(i=0; i<n; i++)
	{
	  temp=fabs(g[i])*MAX(fabs(p[i]),1.0)/den;
	  if(temp>test)test=temp;
	}
      if(test<gtol) return;
      for(i=0; i<n; i++) dg[i]=g[i]-dg[i];
      for(i=0; i<n; i++)
	{
	  hdg[i]=0.0;
	  for(j=0; j<n; j++) hdg[i] += hessin[i][j]*dg[j];
	}
      fac=fae=sumdg=sumxi=0.0;
      for(i=0; i<n; i++)
	{
	  fac += dg[i]*xi[i];
	  fae += dg[i]*hdg[i];
	  sumdg += dg[i]*dg[i];
	  sumxi += xi[i]*xi[i];
	}
      if(fac > sqrt(EPS*sumdg*sumxi))
	{
	  fac=1.0/fac;
	  fad=1.0/fae;
	  for(i=0; i<n; i++) dg[i]=fac*xi[i]-fad*hdg[i];
	  for(i=0; i<n; i++)
	    {
	      for(j=i; j<n; j++)
		{
		  hessin[i][j] += fac*xi[i]*xi[j]
		    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
		  hessin[j][i]=hessin[i][j];
		}
	    }
	}
      for(i=0; i<n; i++)
	{
	  xi[i]=0.0;
	  for(j=0; j<n; j++) xi[i] -= hessin[i][j]*g[j];
	}
    }
}

/* ---------------------------------------------------------------------- */

double SEGPdist::
PairLogEvidence(vector<int>& itemIndex,
			 double& lengthScale,
			 double& noiseFreeScale,
			 double& noiseSigma)
{
  // Declarations

  int i, j, index;

  const int nCurrentItems=itemIndex.size();

  double replicateNoise,logEvidence=-numeric_limits<double>::infinity();

  vector<double> yValues=vector<double>(nCurrentItems*nTimePoints);

  // Extract the data points for this current cluster;
  // store the relevant data in yValues in row-major order

  for (i=0; i<nCurrentItems; i++)
  {
    index=i;
    for (j=0; j<nTimePoints; j++)
    {
      yValues[index]=data[itemIndex[i]][j];
      index+=nCurrentItems;
    }
  }

  // Optimise the hyperparameters (length-scale, noise-free-scale, noise-sigma)

  if (noise_mode == 0 )
  {
    OptimiseHyperparameters(yValues, lengthScale, noiseFreeScale, noiseSigma);
    logEvidence = ComputeMaximisedLogEvidence(yValues, lengthScale, noiseFreeScale, noiseSigma);

  }else{

    cout << "Error: noise_mode not recogised" <<endl;
  }
  return logEvidence;
}

/* ---------------------------------------------------------------------- */

double  SEGPdist::
ComputeMaximisedLogEvidence(vector<double> yValues,
			    double& lengthScale,
			    double& noiseFreeScale,
			    double& noiseSigma)
{
  // Declarations
  int blockSize;
  BlockCovarianceMatrix covarFunction;
  double logEvidence;
  
  // Find block size
  blockSize = yValues.size() / nTimePoints;

  // Calculate maximised log-evidence
  covarFunction = SquareExponentialCovarianceFunction(lengthScale,
						      blockSize, 
						      noiseFreeScale);
  covarFunction = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
  logEvidence = ComputeLogEvidence(covarFunction, yValues);

  return logEvidence;
}

/* ---------------------------------------------------------------------- */

void SEGPdist::
OptimiseHyperparameters(const vector<double>& yValues,
			double& lengthScale,
			double& noiseFreeScale,
			double& noiseSigma)
{
  int i;
  int blockSize = yValues.size() / nTimePoints;
  double bestLengthScale;
  double bestLogEv, trialLogEv;
  vector<double> params(3);
  vector<int> fix(3);
  
  // Guess a starting state using a coarse-grained method
  bestLengthScale=2.0;
  params[1]=1.0;
  params[2]=0.5;
  bestLogEv=-numeric_limits<double>::infinity();
  for(i=2; i<=10; i+=2)
    {
      params[0] = static_cast<double>(i);
      trialLogEv = ComputeLogEvidenceFromHyperparameters(yValues,
							 blockSize,
							 params);
      if(trialLogEv > bestLogEv)
	{
	  bestLengthScale = params[0];
	  bestLogEv = trialLogEv;
	}
    }
  params[0]=bestLengthScale;

  // Now do the actual maximisation
  fix[0]=fix[1]=fix[2]=0; // do not fix any params
  double fret=0;
  double gtol=fast_switch ? (1.0e-1) : (1.0e-2); // the convergence tolerance
  DFPMaximise(params,fix,gtol,fret,blockSize,yValues);
  
  // Return the result
  lengthScale=params[0];
  noiseFreeScale=params[1];
  noiseSigma=params[2];
}

/* ---------------------------------------------------------------------- */

void SEGPdist::
ImposeConstraintsOnHyperparameters(vector<double>& params)
{
  params[0] = MAX(params[0], 0.2); // lengthScale
  params[1] = MAX(params[1], 0.2); // noiseFreeScale
  params[2] = MAX(params[2], 0.05); // noiseSigma
  params[2] = MIN(params[2], 1.0); // noiseSigma
}

/* ---------------------------------------------------------------------- */

double SEGPdist::
ComputeLogEvidenceFromHyperparameters(const vector<double>& yValues,
				      const int blockSize,
				      const vector<double>& params)
{
  BlockCovarianceMatrix covarFunc = 
    AddNoiseToCovarianceFunction(SquareExponentialCovarianceFunction(params[0],
								     blockSize,
								     params[1]),
				 params[2]);
  return ComputeLogEvidence(covarFunc, yValues);
}

/* ---------------------------------------------------------------------- */

void SEGPdist::
ComputeGradientsFromHyperparameters(const vector<double>& yValues,
				    const int blockSize,
				    const vector<double>& params,
				    vector<double>& grad)
{
  // Build the covariance function for these hyperparameters
  BlockCovarianceMatrix covarFunc = 
    AddNoiseToCovarianceFunction(SquareExponentialCovarianceFunction(params[0],
								     blockSize,
								     params[1]),
				 params[2]);
  // Useful intermediate values for computing the gradient
  covarFunc.InvertMatrix(); // K = K^-1
  vector<double> alpha = covarFunc.VectorMultiply(yValues); // K^-1 * y
  BlockCovarianceMatrix covarDeriv_ls =
    SquareExponential_lengthDerivative(params[0],
				       blockSize,
				       params[1]);
  BlockCovarianceMatrix covarDeriv_nf = SquareExponentialCovarianceFunction(params[0],
									    blockSize,
									    1);
  // Compute the gradient at this point in hyperspace
  grad[0] = ComputeGradient(covarFunc, covarDeriv_ls, alpha);
  grad[1] = ComputeGradient(covarFunc, covarDeriv_nf, alpha);
  grad[2] = ComputeNoiseGradient(covarFunc, alpha, params[2]);
}

/* ---------------------------------------------------------------------- */

void SEGPdist::
ComputeLogEvidenceAndGradientsFromHyperparameters(const vector<double>& yValues,
						  const int blockSize,
						  const vector<double>& params,
						  double& logEv,
						  vector<double>& grad)
{
  // Build the covariance function for these hyperparameters
  BlockCovarianceMatrix covarFunc = 
    AddNoiseToCovarianceFunction(SquareExponentialCovarianceFunction(params[0],
								     blockSize,
								     params[1]),
				 params[2]);
  // Compute the log-evidence
  logEv = ComputeLogEvidence(covarFunc, yValues);
  
  // Useful intermediate values for computing the gradient
  covarFunc.InvertMatrix(); // K = K^-1
  vector<double> alpha = covarFunc.VectorMultiply(yValues); // K^-1 * y
  BlockCovarianceMatrix covarDeriv_ls =
    SquareExponential_lengthDerivative(params[0],
				       blockSize,
				       params[1]);
  BlockCovarianceMatrix covarDeriv_nf = SquareExponentialCovarianceFunction(params[0],
									    blockSize,
									    1);
  // Compute the gradient at this point in hyperspace
  grad[0] = ComputeGradient(covarFunc, covarDeriv_ls, alpha);
  grad[1] = ComputeGradient(covarFunc, covarDeriv_nf, alpha);
  grad[2] = ComputeNoiseGradient(covarFunc, alpha, params[2]);
}

/* ---------------------------------------------------------------------- */

BlockCovarianceMatrix SEGPdist::
SquareExponentialCovarianceFunction(double lengthScale,
				    int blockSize,
				    double noiseFreeScale)
{
  // Declarations
  int                    i, j;
  double                 covarElement;
  BlockCovarianceMatrix  blockMatrix;

  // Initialise the block matrix
  blockMatrix.nRank     = nTimePoints;
  blockMatrix.blockSize = blockSize;

  // Initialise the covariance function
  blockMatrix.noisyCoeff=vector<double>(nTimePoints, 0.0);
  blockMatrix.noiseFreeCoeff=
    vector<vector<double> >(nTimePoints, vector<double>(nTimePoints, 0.0));

  // Compute each element of the covariance function
  for (i=0; i<nTimePoints; i++)
  {
    for (j=i; j<nTimePoints; j++)
    {
      covarElement = fabs(timePoints[i] - timePoints[j]);
      covarElement = covarElement*covarElement;
      covarElement /= 2*lengthScale*lengthScale;
      covarElement *= -1;
      covarElement  = exp(covarElement);
      covarElement *= noiseFreeScale;
      // and store in 2 elements (covariance is symmetric)
      // this duplicates effort for i==j; not a big deal, computationally :-)
      blockMatrix.noiseFreeCoeff[i][j] = covarElement;
      blockMatrix.noiseFreeCoeff[j][i] = covarElement;
    }
  }
  
  return blockMatrix;
}

/* ----------------------------------------------------------------------
   Compute the partial derivative w.r.t. length-scale of a square
   exponential (SE) covariance function. Note that we've hard-coded
   noise=0 here.
---------------------------------------------------------------------- */

BlockCovarianceMatrix SEGPdist::
SquareExponential_lengthDerivative(double lengthScale,
				   int blockSize,
				   double noiseFreeScale)
{
  // Declarations

  int i, j;
  double covarElement, deltaTime;
  BlockCovarianceMatrix  blockMatrix;

  assert(lengthScale > 0);
  assert(blockSize > 0);

  // Initialise the block matrix

  blockMatrix.nRank = nTimePoints;
  blockMatrix.blockSize = blockSize;

  // Initialise the covariance function
  blockMatrix.noiseFreeCoeff=
    vector<vector<double> >(nTimePoints, vector<double>(nTimePoints, 0.0));
  blockMatrix.noisyCoeff = vector<double>(nTimePoints, 0.0);

  // Compute each element of the covariance function

  for (i=0; i<nTimePoints; i++)
  {
    for (j=i; j<nTimePoints; j++)
    {
      deltaTime     = fabs(timePoints[i] - timePoints[j]);
      covarElement  = deltaTime*deltaTime;
      covarElement /= 2*lengthScale*lengthScale;
      covarElement *= -1;
      covarElement  = exp(covarElement);
      covarElement *= deltaTime*deltaTime;
      covarElement /= lengthScale*lengthScale*lengthScale;
      covarElement *= noiseFreeScale;
      // and store in 2 elements (covariance is symmetric)
      // this duplicates effort for i==j; not a big deal, computationally :-)
      blockMatrix.noiseFreeCoeff[i][j] = covarElement;
      blockMatrix.noiseFreeCoeff[j][i] = covarElement;
    }
  }

  return blockMatrix;
}

/* ---------------------------------------------------------------------- */
