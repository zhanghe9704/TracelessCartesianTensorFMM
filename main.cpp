#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include "FMMframe.hpp"

#include <cstring>
#include <cmath>

#include <fstream>

using std::cout;
using std::endl;
using std::vector;

int Coulomb(double * x,double * y,double * z,double * q, unsigned long int n_ptc, unsigned long int n_calc,double * phi){

    if(n_calc>n_ptc){
        cout<<"Error! in Coulomb"<<endl;
        return 1;
    }

    memset(phi, 0, n_calc*sizeof(double));

    double r = 0;
    for(unsigned long int i=0; i<n_calc; ++i){
        for(unsigned long int j=0; j<n_ptc; ++j){
            if(i!=j){
                r = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2));
                phi[i] += q[j]/r;
            }
        }
    }

    return 0;
}

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

	unsigned long int N=100;

	double * x = new double[N];
	double * y = new double[N];
	double * z = new double[N];
	double * q = new double[N];
	double * phi = new double[N];
	double * phi_check = new double[N];

//	 Generate random particle positions.
//	srand (time(NULL));
//	for(unsigned long int i=0;i<N;++i){
//		x[i] = rand();
//		y[i] = rand();
//		z[i] = rand();
//		q[i] = 1;
//	}
//
//	scale(x, N);
//	scale(y, N);
//	scale(z, N);
//
	char filename[30] = "xyz.txt";
//    std::ofstream outfile;
//    outfile.open(filename);
//    for(unsigned long int i=0;i<N;++i){
//        outfile<<x[i]<<' '<<y[i]<<' '<<z[i]<<endl;
//    }
//    outfile.close();
//    return 0;

    std::ifstream infile;
    infile.open(filename);
    for(unsigned long int i=0; i<N; ++i){
        infile>>x[i];
        infile>>y[i];
        infile>>z[i];
        q[i] = 1;
    }
    infile.close();

    cout<<x[0]<<' '<<y[0]<<' '<<z[0]<<endl;
    cout<<x[N-1]<<' '<<y[N-1]<<' '<<z[N-1]<<endl;



    Coulomb(x,y,z,q,N,N,phi_check);

    int n_ptc_box = 10;
    int n_rank = 6;

    fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);

    char filename2[30] = "checkphi.txt";
    std::ofstream output;
    output.open(filename2);
    for(unsigned long int i=0;i<N;++i){
        output<<phi[i]<<' '<<phi_check[i]<<' '<<phi[i]-phi_check[i]<<endl;
    }
    output.close();

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] q;
    delete[] phi;
    delete[] phi_check;
	system("pause");
	return 0;
}
