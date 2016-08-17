#ifndef GRAD_H
#define GRAD_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <math.h>

#include "stringtools.h"
#include "pTable.h"

class Gradient {
	private:
   int natoms;
   int N3;
   string* anames;
   int* anumbers;
	double* amasses; 
	public: 
   	double* grad; 
		void grad_read(string gradfile,double* grad);	
};

#endif
