#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>

#include "FMMframe.hpp"

#include <cstring>
#include <cmath>
#include <random>
#include <chrono>

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

double error(double * phi, double * phi_fmm, unsigned long int n){
    double sum_phi = 0;
    double sum_dphi = 0;
    for(unsigned long int i=0; i<n; ++i){
        sum_phi += pow(phi[i],2);
        sum_dphi += pow(phi[i]-phi_fmm[i],2);
    }
    return sqrt(sum_dphi/sum_phi);
}


int main(){

	unsigned long int N=256000;
	unsigned long int N_calc = 1000;
	int n_ptc_box = 1000;
    int n_rank = 6;

	double * x = new double[N];
	double * y = new double[N];
	double * z = new double[N];
	double * q = new double[N];
	double * phi = new double[N];
	double * phi_check = new double[N];

////	 Generate random particle positions.
	srand (time(NULL));
	for(unsigned long int i=0;i<N;++i){
		x[i] = rand();
		y[i] = rand();
		z[i] = rand();
		q[i] = 1;
	}

////  normal distribution.
//    // obtain a seed from the timer
//    std::default_random_engine generator;
//    generator.seed(time(NULL));
//    std::normal_distribution<double> distribution(0.0,1.0);
//    for(unsigned long int i=0;i<N;++i){
//        x[i] = distribution(generator);
//        y[i] = distribution(generator);
//        z[i] = distribution(generator);
//        q[i] = 1;
//    }

	scale(x, N);
	scale(y, N);
	scale(z, N);
	cout<<N<<" particles initialized!"<<endl;
//
//	char filename[30] = "xyz.txt";
////    std::ofstream outfile;
////    outfile.open(filename);
////    for(unsigned long int i=0;i<N;++i){
////        outfile<<x[i]<<' '<<y[i]<<' '<<z[i]<<endl;
////    }
////    outfile.close();
////    return 0;
//
//    std::ifstream infile;
//    infile.open(filename);
//    for(unsigned long int i=0; i<N; ++i){
//        infile>>x[i];
//        infile>>y[i];
//        infile>>z[i];
//        q[i] = 1;
//    }
//    infile.close();

    std::chrono::steady_clock::time_point start, end;

    cout<<"Calculate the first "<<N_calc<<" particles by the pairwise Coulomb formula!"<<endl;
    start = std::chrono::steady_clock::now();
    Coulomb(x,y,z,q,N,N_calc,phi_check);
    end = std::chrono::steady_clock::now();
    auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout<<"Coulomb formula: "<<t1<<' '<<N/N_calc*t1<<endl;

    cout<<"start FMM"<<endl;
    start = std::chrono::steady_clock::now();
    fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
    end = std::chrono::steady_clock::now();
    auto t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout<<"FMM: "<<t3<<endl;
    cout<<"end FMM"<<endl;
    cout<<"FMM/Coulomb: "<<t3/(N/N_calc*t1)<<endl;

//    cout<<"Calculate "<<N<<" particles by the pairwise Coulomb formula!"<<endl;
//    start = std::chrono::steady_clock::now();
//    Coulomb(x,y,z,q,N,N,phi_check);
//    end = std::chrono::steady_clock::now();
//    auto t2 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//    cout<<"Coulomb formula: "<<t2<<endl;
//
//
//    cout<<"FMM/Coulomb: "<<t3/t2<<endl;

//    char filename2[30] = "checkphi.txt";
//    std::ofstream output;
//    output.open(filename2);
//    for(unsigned long int i=0;i<N_calc;++i){
//        output<<phi[i]<<' '<<phi_check[i]<<' '<<phi[i]-phi_check[i]<<endl;
//    }
//    output.close();

    cout<<"Relative error: "<<error(phi,phi_check,N_calc)<<endl;

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] q;
    delete[] phi;
    delete[] phi_check;
	system("pause");
	return 0;
}
