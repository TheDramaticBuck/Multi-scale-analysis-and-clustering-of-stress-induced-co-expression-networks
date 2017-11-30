/* ----------------------------------------------------------------------
   
	Compute Gaussian Process inspired distance matrix


 	For specific details on the underlying methods and the datasets see:
 
        https://arxiv.org/abs/1703.02872 


	Dr. Nuno R. Nene, University of Cambridge, 2017

 	Email: nunonene@gmail.com
   
------------------------------------------------------------------------- */


/* 

Notes: This code was used to calculate a dissimilarity matrix between 

expression profiles which had no replicates. If replicated experiments

are included, adjustments need to be made.

*/


#if defined(macintosh) || defined(MACOS) || defined(_MACOS)  || defined(__APPLE__)  || defined(_MAC)   || defined(MAC)  || defined(mac)  || defined(MACINTOSH)
#define __MAC__
#endif


#if defined(_WIN32)  || defined(WIN32)  || defined(__WIN32__) || defined(_WIN64)  || defined(WIN64) || defined(__WIN64__) || defined(_WINDOWS) || defined(__WINDOWS__)
#define __WIN__
#endif


#if defined(linux)  || defined(__linux__)  || defined(_linux) || defined(LINUX)  || defined(_LINUX) || defined(_UNIX) || defined(__UNIX__) || defined(__gnu_linux__) || defined(__unix__) || defined(UNIX) || defined(unix) || defined(sparc)
#define __lin__
#endif





#ifdef __WIN__
#include <time.h>
#include <sys/timeb.h> 
#include <process.h>
#include <windows.h>
#endif


#include <limits>

#include <stdio.h>

#include "matrix.h"

#include "header.h"

#include "GPdist.h"

#include "SEGPdist.h"

#include "BlockCovarianceMatrix.h"


extern "C" {

#include <math.h>
#include "mex.h"


}


#ifdef __WIN__
int gettimeofday (struct timeval *tp, void *tz)
{
	struct _timeb timebuffer;
	_ftime (&timebuffer);
	tp->tv_sec = timebuffer.time;
	tp->tv_usec = timebuffer.millitm * 1000;
	return 0;
}
#endif


#ifdef __lin__
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h> 
#endif

#ifdef __MAC__
#include <sys/time.h>
#endif


extern "C" { 

  void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
     
  
	// Declarations

	int i, j, counter1=0;
    
	vector<double> GPmat; // vector for storage of dissimilarity matrix in row-form
    
	SEGPdist* segp=NULL; // Class where most of the routines concerning a Gaussian Process with a Squared-exponential are defined
    
	GPdist GP; // Class linking to SEGPdist
       
    	vector<double> noise, timepoints_copy;

    	vector<vector<double>> data_mat; // time-series data
    
    
    
    	if (nrhs>0){
            
		//Inputs:

        	double *data = (double*) mxGetPr(prhs[0]);             // 0. double* inputData 

        	double *timepoints = (double*) mxGetPr(prhs[1]);       // 1. double* timepoints

		int nGenes=(int)mxGetScalar(prhs[2]);                  // 2. int* nGenes

        	int nTimepoints=(int)mxGetScalar(prhs[3]);             // 3. int* nTimepoints
             


     		// Copy data 

     		for (i=0; i<nGenes; i++){

        	 data_mat.push_back(vector<double>(nTimepoints, 0));
          
        	 	for (j=0; j<nTimepoints; j++){
           		data_mat[i][j] = data[counter1++];
         		}
     		}

        	//Read in the timepoints

     		for (i=0; i<nTimepoints; i++){

          		timepoints_copy.push_back(timepoints[i]);
        	}


        	// Instantiate the required class
       	 
		segp = new SEGPdist(data_mat);

		noise.push_back(0.0);
           	 
		segp->ReadInNoise(noise);

        	segp->SetNoiseMode(0);

        	segp->ReadInTimePoints(timepoints_copy);

		segp->SetReps(1);


		// Compute GP dissimilarity matrix
   
        	GPmat = GP.ComputeGPdist(*segp); // output is the log of the ratio between joint and the product of separate probabilities; the distance matrix is 							    computed in Matlab


    	 	// Cleanup

   	  	delete segp;
		
		data_mat.clear();

		noise.clear();

		timepoints_copy.clear();

     
    	 	//Output (into Matlab)
     
    	 	plhs[0] = mxCreateDoubleMatrix(nGenes*(nGenes-1)/2, 1, mxREAL);

			
     		double * output_GPdist = (double*) mxGetPr(plhs[0]);
     
     
     		for (int j = 0; j <(nGenes*(nGenes-1)/2) ; j++) {
				 
			output_GPdist[j]=(double) GPmat[j];
     		}

		GPmat.clear();
     
    	}else{
	 	cout<<"Number of arguments incorrect"<<endl;
		return;
    	}
    
  }

    
}
