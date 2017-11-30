/* ----------------------------------------------------------------------
	
	Adapted by:

 	Dr. Nuno R. Nene, University of Cambridge, 2017

 	Email: nunonene@gmail.com


 	For specific details on the underlying methods and the datasets see:
 
        https://arxiv.org/abs/1703.02872 


	Based on/Adapted from C++ files provided with:
   
	http://www.bioconductor.org/packages/release/bioc/html/BHC.html

   
------------------------------------------------------------------------- */


#include <numeric>



#include "BlockCovarianceMatrix.h"

extern "C" {
#include <math.h>
#include "mex.h"
}

/* ---------------------------------------------------------------------- */

BlockCovarianceMatrix::BlockCovarianceMatrix() {}

/* ---------------------------------------------------------------------- */

BlockCovarianceMatrix::
BlockCovarianceMatrix(const double noisyValue,
		      const double noiseFreeValue,
		      const int inputBlockSize)
{
  nRank = 1;
  blockSize = inputBlockSize;
  noisyCoeff.push_back(noisyValue);
  noiseFreeCoeff.push_back(vector<double>(1, noiseFreeValue));
}

/* ---------------------------------------------------------------------- */


BlockCovarianceMatrix BlockCovarianceMatrix::
Build_E_SubMatrix(const BlockCovarianceMatrix& inputMatrix)
{
  // The object we are going to construct
  BlockCovarianceMatrix bcm;

  // Declarations
  double noisyScaling, newNoiseFree, temp2_newNoiseFree;
  BlockCovarianceMatrix subMatrix_A;

  // This won't work with a rank-1 input matrix
  assert(inputMatrix.nRank > 1);

  // Assign some obvious values
  bcm.nRank = inputMatrix.nRank - 1;
  bcm.blockSize = inputMatrix.blockSize;

  // Find the sub-matrix "A" and invert it
  subMatrix_A = BlockCovarianceMatrix(inputMatrix.noisyCoeff[0],
				      inputMatrix.noiseFreeCoeff[0][0],
				      bcm.blockSize);
  subMatrix_A.InvertRankOneMatrix();

  // Initialise the arrays
  bcm.noisyCoeff = vector<double>(bcm.nRank);
  bcm.noiseFreeCoeff=vector<vector<double> >(bcm.nRank, vector<double>(bcm.nRank));

  // Assign the sub-matrix values
  const double temp_newNoiseFree = (subMatrix_A.noisyCoeff[0] + bcm.blockSize)*subMatrix_A.noiseFreeCoeff[0][0]*bcm.blockSize;
  for (int i=0; i<bcm.nRank; i++)
  {
    temp2_newNoiseFree=temp_newNoiseFree*inputMatrix.noiseFreeCoeff[i+1][0]; // B
    for (int j=0; j<bcm.nRank; j++)
    {
      // find adjusting noiseFree term
      newNoiseFree = temp2_newNoiseFree * inputMatrix.noiseFreeCoeff[0][j+1]; // C

      // make adjustments
      bcm.noiseFreeCoeff[i][j] = inputMatrix.noiseFreeCoeff[i+1][j+1] - newNoiseFree;
    }
    // find and adjust the noisy term for this row
    noisyScaling  = inputMatrix.noiseFreeCoeff[i+1][i+1] / bcm.noiseFreeCoeff[i][i];
    bcm.noisyCoeff[i] = inputMatrix.noisyCoeff[i+1] * noisyScaling;
  }

  return bcm;
}


/* ---------------------------------------------------------------------- */

void BlockCovarianceMatrix::InvertMatrix()
{
  if(nRank==1)
  {
    InvertRankOneMatrix();
  }
  else
  {
    InvertBlockMatrix();
  }
}

/* ---------------------------------------------------------------------- */

void BlockCovarianceMatrix::InvertBlockMatrix()
{
  // Declarations
  int i, j;
  BlockCovarianceMatrix subMatrix_A, subMatrix_E;
  vector<double>        noiseFree_B, offDiagonal;
  double                diagonal_noisy, diagonal_noiseFree, factor_AB;

  // Use the block-matrix structure to our advantage;
  // find the sub-matrices A and E
  subMatrix_A = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize);
  subMatrix_E = Build_E_SubMatrix(*this);

  // Invert A and E
  subMatrix_A.InvertRankOneMatrix();
  subMatrix_E.InvertMatrix();

  // Find the sub-matrices B and C
  // (symmetry => only need one of these)
  noiseFree_B = noiseFreeCoeff[0];
  // remove the first element as this is part of A, not B
  noiseFree_B.erase(noiseFree_B.begin(), noiseFree_B.begin()+1);

  // Find the factor coming from (A^-1)*B
  factor_AB  = subMatrix_A.noiseFreeCoeff[0][0];
  factor_AB *= blockSize + subMatrix_A.noisyCoeff[0];

  // Construct the overall inverse matrix;
  // copy E^-1 into the D slots
  for (i=1; i<nRank; i++)
  {
    noisyCoeff[i] = subMatrix_E.noisyCoeff[i-1];
    for (j=1; j<nRank; j++)
    {
      noiseFreeCoeff[i][j] = subMatrix_E.noiseFreeCoeff[i-1][j-1];
    }
  }
  // construct the remaining off-diagonal elements, using B * E^-1
  offDiagonal = subMatrix_E.BlockMultiply(noiseFree_B);
  for (i=1; i<nRank; i++)
  {
    // using the fact that our matrix is symmetric
    noiseFreeCoeff[0][i] = noiseFreeCoeff[i][0] = -factor_AB * offDiagonal[i-1];
  }
  // construct the final diagonal element
  diagonal_noiseFree = 0;
  for (i=0; i<(nRank-1); i++) // need to find B * E^-1 * C here...
  {
    // using the fact that C_transpose = B
    diagonal_noiseFree += blockSize * offDiagonal[i] * noiseFree_B[i];
  }
  diagonal_noiseFree *= factor_AB * factor_AB;
  // also need to add on A^-1
  diagonal_noiseFree += subMatrix_A.noiseFreeCoeff[0][0];
  // hence compute the new noisy term
  diagonal_noisy = subMatrix_A.noisyCoeff[0] * subMatrix_A.noiseFreeCoeff[0][0];
  diagonal_noisy /= diagonal_noiseFree;
  // and store the values in this object
  noiseFreeCoeff[0][0] = diagonal_noiseFree;
  noisyCoeff[0] = diagonal_noisy;
}

/* ---------------------------------------------------------------------- */
void BlockCovarianceMatrix::InvertRankOneMatrix()
{
  // Declarations
  double newNoisy, newNoiseFree;

  // Compute the new noisy value (for the diagonal elements)
  newNoisy  = -noisyCoeff[0] - blockSize;

  // Compute the new noise-free value
  newNoiseFree  = -1.0 /
    (noiseFreeCoeff[0][0] * (noisyCoeff[0] * (noisyCoeff[0] + blockSize)));

  // Update the matrix values, so this matrix object is now inverted
  noisyCoeff[0] = newNoisy;
  noiseFreeCoeff[0][0] = newNoiseFree;
}


/* ---------------------------------------------------------------------- */

double BlockCovarianceMatrix::ComputeMatrixDeterminant() const
{
  // Declarations
  double logDeterminant;
  BlockCovarianceMatrix subMatrix_A, subMatrix_E;//use the maths notation here

  // Recursion to find the log-det
  if (nRank==1)
  {
    logDeterminant = ComputeRankOneMatrixDeterminant();
  }
  else
  {
    // extract rank n-1 and rank-1 sub-matrices
    // remember that we want "E"
    subMatrix_A = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize);
    subMatrix_E = Build_E_SubMatrix(*this);

    // compute the contributions from the sub-matrices
    logDeterminant  = subMatrix_A.ComputeRankOneMatrixDeterminant();
    logDeterminant += subMatrix_E.ComputeMatrixDeterminant();
  }

  return logDeterminant;
}

/* ---------------------------------------------------------------------- */

double BlockCovarianceMatrix::ComputeRankOneMatrixDeterminant() const
{
  // Declarations
  double logDeterminant;

  assert(nRank==1);

  if (noiseFreeCoeff[0][0] < 0.0 || noisyCoeff[0] < 0.0)
  {
    cout << "problem with log Determinant: nan" << endl;
  }

  //Compute the log-det
  logDeterminant  = log(noiseFreeCoeff[0][0]) * blockSize;
  logDeterminant += log(noisyCoeff[0]) * (blockSize - 1);
  logDeterminant += log(noisyCoeff[0] + blockSize);

  return logDeterminant;
}

/* ---------------------------------------------------------------------- */

vector<double> BlockCovarianceMatrix::
BlockMultiply(const vector<double>& inputVector) const
{
  // Declarations
  int i;
  double currentElement;
  const int input_size=inputVector.size();
  vector<double> outputVector=vector<double>(input_size);

  // Compute the output vector elements
  for (i=0; i<input_size; i++)
  {
    // Inner product: <noiseFreeCoeff row, inputVector>
    currentElement = std::inner_product(inputVector.begin(),
					inputVector.end(),
					noiseFreeCoeff[i].begin(),
					0.0);

    // Since it's a block matrix, we use each matrix term blockSize-many times
    currentElement *= blockSize;

    // Add the (diagonal) noise term
    currentElement += noisyCoeff[i] * noiseFreeCoeff[i][i] * inputVector[i];

    // Store the result
    outputVector[i] = currentElement;
  }

  return outputVector;
}

/* ---------------------------------------------------------------------- */

double BlockCovarianceMatrix::
ComputeLogLikelihoodProduct(const vector<double>& data) const
{
  // Compute K * y
  vector<double> temp=VectorMultiply(data);
  
  // Compute and return y^t * (K * y)
  return std::inner_product(temp.begin(), temp.end(), data.begin(), 0.0);
}

/* ---------------------------------------------------------------------- */

double BlockCovarianceMatrix::GetElement(const int i, const int j) const
{
  // Declarations
  int block_i, block_j;
  double matrixElement;

  // Find which block the element belongs to
  block_i = i / blockSize;
  block_j = j / blockSize;

  // Find the matrix element
  matrixElement = noiseFreeCoeff[block_i][block_j];

  // If it's a diagonal element we make an adjustment
  if (i==j) matrixElement *= 1.0 + noisyCoeff[block_i];

  return matrixElement;
}

/* ---------------------------------------------------------------------- */

vector<double> BlockCovarianceMatrix::
VectorMultiply(const vector<double>& inputVector) const
{
  // Declarations
  vector<double> outputVector=vector<double>(inputVector.size());
  vector<double> row;
  vector<double>::const_iterator inIt, inIt_e;
  int i,j, out_id,out_id_e;
  double dot, noisyFactor, temp;

  out_id=0;
  for(i=0; i<nRank; i++)
    {
      row = noiseFreeCoeff[i];
      noisyFactor = noisyCoeff[i]*row[i];
      inIt = inputVector.begin();
      dot=0;
      
      // Iterate over each block in this row
      for(j=0; j<nRank; j++)
	{
	  // Apply the coeff to blockSize-many elements from the input vector
	  inIt_e = inIt+blockSize;
	  temp = std::accumulate(inIt, inIt_e, 0.0);
	  dot += temp*row[j];

	  // Move our iterator to the next block of our vector
	  inIt = inIt_e;
	}
      
      // Store the result
      out_id_e = out_id + blockSize;
      for(; out_id<out_id_e; out_id++)
	{
	  outputVector[out_id] = dot + inputVector[out_id]*noisyFactor;
	}
    }
  
  return outputVector;
}

/* ---------------------------------------------------------------------- */

