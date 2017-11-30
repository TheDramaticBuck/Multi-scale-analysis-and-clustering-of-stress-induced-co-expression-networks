/* ----------------------------------------------------------------------
 C++ files provided with:
   
	http://www.bioconductor.org/packages/release/bioc/html/BHC.html
------------------------------------------------------------------------- */

#ifndef HEADER_H
#define HEADER_H

//#define NDEBUG

// Handy macros
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

extern bool fast_switch; // declared in header.cpp

static const double dirichletProcessParameter=0.001;

// Standard includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <climits>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <numeric>
#include <ctime>


#include "gammaln.h"

using namespace std;

#endif // HEADER_H

