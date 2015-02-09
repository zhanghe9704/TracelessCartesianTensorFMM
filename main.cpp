#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include "FMMframe.hpp"

using std::cout;
using std::endl;
using std::vector;


int scale(double *x, unsigned long int n){
	double max_x = x[0];
	double min_x = x[0];
	for(unsigned long int i=0; i<n; ++i){
		if(max_x<x[i])  max_x = x[i];
		if(min_x>x[i])  min_x = x[i];
	}

	double avg = 0.5*(max_x+min_x);
	double scale = 2/(max_x-min_x);
	for (unsigned long int i=0; i<n; ++i){
		x[i] = scale*(x[i]-avg);
	}

	return 0;
}

int main(){

	unsigned long int N=40;

	double * x = new double[N];
	double * y = new double[N];
	double * z = new double[N];
	double * q = new double[N];
	double * phi = new double[N];

//	 Generate random particle positions.
	srand (time(NULL));
	for(unsigned long int i=0;i<N;++i){
		x[i] = rand();
		y[i] = rand();
		z[i] = rand();
		q[i] = 1;
	}

	scale(x, N);
	scale(y, N);
	scale(z, N);



    fmm(x,  y,  z,  q, N, 4, 4, phi);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] q;
    delete[] phi;
	system("pause");
	return 0;
}
