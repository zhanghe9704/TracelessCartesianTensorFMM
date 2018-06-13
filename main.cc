#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "fmm_frame.h"

using std::cout;
using std::endl;
using std::vector;

int Coulomb(double * x,double * y,double * z,double * q, unsigned long int n_ptc, unsigned long int n_calc,double * phi){

    if(n_calc>n_ptc){
        cout<<"Error! in Coulomb"<<endl;
        return 1;
    }

    memset(phi, 0, n_calc*sizeof(double));

    for(unsigned long int i=0; i<n_calc; ++i){
        for(unsigned long int j=0; j<n_ptc; ++j){
            if(i!=j){
                double dx = x[i]-x[j];
                double dy = y[i]-y[j];
                double dz = z[i]-z[j];
                double r = sqrt(dx*dx+dy*dy+dz*dz);
                phi[i] += q[j]/r;
            }
        }
    }

    return 0;
}

int Coulomb(double * x,double * y,double * z,double * q, unsigned long int n_ptc, unsigned long int n_calc,double * Ex, double * Ey, double * Ez){

    if(n_calc>n_ptc){
        cout<<"Error! in Coulomb"<<endl;
        return 1;
    }

    memset(Ex, 0, n_calc*sizeof(double));
    memset(Ey, 0, n_calc*sizeof(double));
    memset(Ez, 0, n_calc*sizeof(double));

    for(unsigned long int i=0; i<n_calc; ++i){
        for(unsigned long int j=0; j<n_ptc; ++j){
            if(i!=j){
                double dx = x[i]-x[j];
                double dy = y[i]-y[j];
                double dz = z[i]-z[j];
                double r = dx*dx+dy*dy+dz*dz;
                r = 1/sqrt(r);
                r = r*r*r;
                Ex[i] += q[j]*dx*r;
                Ey[i] += q[j]*dy*r;
                Ez[i] += q[j]*dz*r;
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

double error(double * phi, double * phi_fmm, unsigned long int n){
    double sum_phi = 0;
    double sum_dphi = 0;
    for(unsigned long int i=0; i<n; ++i){
        sum_phi += pow(phi[i],2);
        sum_dphi += pow(phi[i]-phi_fmm[i],2);
    }
    return sqrt(sum_dphi/sum_phi);
}

double error(double * Ex, double * Ey, double * Ez, double * Ex_fmm, double * Ey_fmm, double * Ez_fmm, unsigned long int n){
    double sum_phi = 0;
    double sum_dphi = 0;
    for(unsigned long int i=0; i<n; ++i){
        sum_phi += Ex[i]*Ex[i]+Ey[i]*Ey[i]+Ez[i]*Ez[i];
        double dex = Ex[i]-Ex_fmm[i];
        double dey = Ey[i]-Ey_fmm[i];
        double dez = Ez[i]-Ez_fmm[i];
        sum_dphi += dex*dex+dey*dey+dez*dez;
    }
    return sqrt(sum_dphi/sum_phi);
}


int main(){

//	unsigned long int N=131072;  //
//    unsigned long int N=1048576;  //1M 8^6*4
//    unsigned long int N=2097152;    //2M 8^7
//    unsigned long int N=4194304;    //4M 8^7*2
//    unsigned long int N=8388608;    //8M 8^7*4
//    unsigned long int N=16777216;  //16M 8^8
    unsigned long int N=100000;
	unsigned long int N_calc = 1000;
	int n_ptc_box = 64;
    int n_rank = 4;

	double * x = new double[N];
	double * y = new double[N];
	double * z = new double[N];
	double * q = new double[N];
	double * phi = new double[N];
	double * phi_check = new double[N];
	double * Ex = new double[N];
	double * Ey = new double[N];
	double * Ez = new double[N];
	double * Ex_check = new double[N];
	double * Ey_check = new double[N];
	double * Ez_check = new double[N];

////	 Generate random particle positions.
// Uniform distribution.
//	srand (time(NULL));
//	for(unsigned long int i=0;i<N;++i){
//		x[i] = rand();
//		y[i] = rand();
//		z[i] = rand();
//		q[i] = 1;
//	}
//
//		scale(x, N);
//        scale(y, N);
//        scale(z, N);


//  normal distribution.
//////     obtain a seed from the timer
    std::default_random_engine generator;
    generator.seed(time(NULL));
//    generator.seed(1);
    std::normal_distribution<double> distribution(0.0,1.0);
    for(unsigned long int i=0;i<N;++i){
        x[i] = distribution(generator);
        y[i] = distribution(generator);
        z[i] = distribution(generator);
        q[i] = distribution(generator);
    }
//
	scale(x, N);
	scale(y, N);
	scale(z, N);
	scale(q, N);
	cout<<N<<" particles initialized!"<<endl;

    fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi, Ex, Ey, Ez);
	Coulomb(x,y,z,q,N,N_calc, Ex_check, Ey_check, Ez_check);
	cout<<error(Ex_check,Ey_check, Ez_check,Ex, Ey, Ez,N_calc)<<endl;
	Coulomb(x,y,z,q,N,N_calc, phi_check);
	cout<<error(phi_check,phi,N_calc)<<endl;

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] q;
    delete[] phi;
    delete[] phi_check;
    delete[] Ex;
    delete[] Ey;
    delete[] Ez;
    delete[] Ex_check;
    delete[] Ey_check;
    delete[] Ez_check;
	system("pause");
	return 0;
}
