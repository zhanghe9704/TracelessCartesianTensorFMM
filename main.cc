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


double error2(double * Ex, double * Ey, double * Ez, double * Ex_fmm, double * Ey_fmm, double * Ez_fmm, unsigned long int n){
    double sum_phi = 0;
    double sum_dphi = 0;
    for(unsigned long int i=0; i<n; ++i){
        sum_phi += Ex[i]*Ex[i]+Ey[i]*Ey[i]+Ez[i]*Ez[i];
//        double dex = Ex[i]-Ex_fmm[i];
//        double dey = Ey[i]-Ey_fmm[i];
//        double dez = Ez[i]-Ez_fmm[i];
        double de2 = Ex[i]*Ex[i]+Ey[i]*Ey[i]+Ez[i]*Ez[i]-(Ex_fmm[i]*Ex_fmm[i]+Ey_fmm[i]*Ey_fmm[i]+Ez_fmm[i]*Ez_fmm[i]);
        sum_dphi += fabs(de2);
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
    unsigned long int N=16384;
	unsigned long int N_calc = 1000;
	int n_ptc_box = 64;
    int n_rank = 6;

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
//
////	 Generate random particle positions.
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
//
//
//	for(int i=0; i<N-1000; ++i){
//        if(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]>1) q[i]=0;
//	}
//
//	for(int i=0; i<1000; ++i){
//        x[N-1000+i] = 0;
//        y[N-1000+i] = 0;
//        z[N-1000+i] = (i+1)*0.005;
//        q[N-1000+i] = 0;
//	}
//
//	fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
//
//	std::ofstream outfile;
//	outfile.open("sphere_potential.txt");
//	for(int i=0; i<1000; ++i){
//        outfile<<x[N-1000+i]<<' '<<y[N-1000+i]<<' '<<z[N-1000+i]<<' '<<phi[N-1000+i]<<endl;
//	}
//	outfile.close();



//  normal distribution.
//////     obtain a seed from the timer
    std::default_random_engine generator;
//    generator.seed(time(NULL));
    generator.seed(1);
    std::normal_distribution<double> distribution(0.0,1.0);
    for(unsigned long int i=0;i<N;++i){
        x[i] = distribution(generator);
        y[i] = distribution(generator);
        z[i] = distribution(generator);
        q[i] = 1;
    }
//
	scale(x, N);
	scale(y, N);
	scale(z, N);
	cout<<N<<" particles initialized!"<<endl;

//    fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//	Coulomb(x,y,z,q,N,N_calc, Ex_check, Ey_check, Ez_check);
//	cout<<error(Ex_check,Ey_check, Ez_check,Ex, Ey, Ez,N_calc)<<endl;


//	fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
//	Coulomb(x,y,z,q,N,N_calc, phi_check);
//	cout<<error(phi_check,phi,N_calc)<<endl;

    fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi, Ex, Ey, Ez);
	Coulomb(x,y,z,q,N,N_calc, Ex_check, Ey_check, Ez_check);
	cout<<error(Ex_check,Ey_check, Ez_check,Ex, Ey, Ez,N_calc)<<endl;
	Coulomb(x,y,z,q,N,N_calc, phi_check);
	cout<<error(phi_check,phi,N_calc)<<endl;
    return 0;
//	Coulomb(x,y,z,q,N,N_calc,Ex_check,Ey_check,Ez_check);


//
////
//	char filename_xyz_out[30] = "xyz_100.txt";
//    std::ofstream outfile;
//    outfile.open(filename_xyz_out);
//    for(unsigned long int i=0;i<N;++i){
//        outfile<<x[i]<<' '<<y[i]<<' '<<z[i]<<endl;
//    }
//    outfile.close();

//    return 0;


//////load particles from file
//    char filename_xyz[30] = "xyz.txt";
//    std::ifstream infile;
//    infile.open(filename_xyz);
//    cout<<"start to read data ..."<<endl;
//    for(unsigned long int i=0; i<N; ++i){
//        infile>>x[i];
//        infile>>y[i];
//        infile>>z[i];
//        q[i] = 1;
//    }
//    infile.close();
//    cout<<N<<" particles initialized!"<<endl;



//////============================================
//////Test speed for different configuration
//////============================================
//
//    std::ofstream outfile;
//    outfile.open("0527_speed_Gaussian_1048576.txt");
//    outfile<<"P s t_E t_phi err_E err_phi t_dir_E t_dir_phi t_rate_E t_rate_phi"<<endl;
//    ////  normal distribution.
//////////     obtain a seed from the timer
//    std::default_random_engine generator;
//    generator.seed(time(NULL));
//    std::normal_distribution<double> distribution(0.0,1.0);
//    for(unsigned long int i=0;i<N;++i){
//        x[i] = distribution(generator);
//        y[i] = distribution(generator);
//        z[i] = distribution(generator);
//        q[i] = 1;
//    }
////
//	scale(x, N);
//	scale(y, N);
//	scale(z, N);
//	cout<<N<<" particles initialized!"<<endl;
//
//	std::chrono::steady_clock::time_point start, end;
//    cout<<"Calculate the first "<<N_calc<<" particles by the pairwise Coulomb formula!"<<endl;
//    start = std::chrono::steady_clock::now();
//    Coulomb(x,y,z,q,N,N_calc,Ex_check,Ey_check,Ez_check);
//    end = std::chrono::steady_clock::now();
//    auto t1_field = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//    cout<<"Coulomb formula: "<<t1_field<<' '<<N/N_calc*t1_field<<endl;
//
//    start = std::chrono::steady_clock::now();
//    Coulomb(x,y,z,q,N,N_calc,phi_check);
//    end = std::chrono::steady_clock::now();
//    auto t1_phi = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//    cout<<"Coulomb formula: "<<t1_phi<<' '<<N/N_calc*t1_phi<<endl;
//
//    int ptc_box[7] = {16,32,64,128,256,512,1024};
//    int order_fmm[9] = {2,3,4,5,6,7,8,9,10};
//
//    int Ns = 7;
//    int Np = 9;
////    int ptc_box[3] = {32,64,128};
////    int order_fmm[2] = {2,3};
//
//    for(int i=0; i<Np; ++i){
//            n_rank = order_fmm[i];
//        for(int j=0; j<Ns; ++j){
//            n_ptc_box = ptc_box[j];
//            start = std::chrono::steady_clock::now();
//            fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//            end = std::chrono::steady_clock::now();
//            auto t3_field = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3_field<<endl;
//            cout<<"FMM/Coulomb: "<<t3_field/(N/N_calc*t1_field)<<endl;
//            double tmp_field = error(Ex_check,Ey_check,Ez_check,Ex,Ey,Ez,N_calc);
//
//            start = std::chrono::steady_clock::now();
//            fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
//            end = std::chrono::steady_clock::now();
//            auto t3_phi = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3_phi<<endl;
//            cout<<"FMM/Coulomb: "<<t3_phi/(N/N_calc*t1_phi)<<endl;
//            double tmp_phi = error(phi_check,phi,N_calc);
//
//
//            outfile<<n_rank<<' '<<n_ptc_box<<' '<<t3_field<<' '<<t3_phi<<' '<<tmp_field<<' '<<tmp_phi<<' '<<t1_field<<' '<<t1_phi<<' '
//            <<t3_field/(N/N_calc*t1_field)<<' '<<t3_phi/(N/N_calc*t1_phi)<<endl;
//        }
//    }
//
//    outfile.close();


////===================================
////Find good s
////===================================

    char filename[50] = "0610_find_s_1e5_rank_2.txt";
    std::ofstream outfile;
    outfile.open(filename);
    outfile<<"N P S T_fmm_phi T_fmm_E"<<endl;

    N = 1e5;
    n_rank = 2;
    int s[7] = {8,16,32,64,128,256,512};
    int ns = 7;

    for(int i=0; i<ns; ++i){

        n_ptc_box = s[i];

    //////            	 Generate random particle positions.
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

        std::chrono::steady_clock::time_point start, end;
        start = std::chrono::steady_clock::now();
        fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
        end = std::chrono::steady_clock::now();
        auto t_phi = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        start = std::chrono::steady_clock::now();
        fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
        end = std::chrono::steady_clock::now();
        auto t_E = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        outfile<<"uniform "<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t_phi<<' '<<t_E<<endl;
        cout<<i<<endl;

    }

        for(int i=0; i<ns; ++i){

        n_ptc_box = s[i];

////         normal distribution.
        ////     obtain a seed from the timer
        std::default_random_engine generator;
        generator.seed(time(NULL));
        std::normal_distribution<double> distribution(0.0,1.0);
        for(unsigned long int i=0;i<N;++i){
            x[i] = distribution(generator);
            y[i] = distribution(generator);
            z[i] = distribution(generator);
            q[i] = 1;
        }

        scale(x, N);
        scale(y, N);
        scale(z, N);

        std::chrono::steady_clock::time_point start, end;
        start = std::chrono::steady_clock::now();
        fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
        end = std::chrono::steady_clock::now();
        auto t_phi = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        start = std::chrono::steady_clock::now();
        fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
        end = std::chrono::steady_clock::now();
        auto t_E = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        outfile<<"uniform "<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t_phi<<' '<<t_E<<endl;
        cout<<i<<endl;

    }

    outfile.close();


//
//////===========================
////// Data generation 06/09/2015
//////===========================
//
//    char filename[50] = "0609_data.txt";
//    std::ofstream outfile;
//    outfile.open(filename);
//    outfile<<"N P S Error_phi T_fmm_phi T_dir_phi Error_E T_fmm_E T_dir_E "<<endl;
//
//    double n_ptc[20] = {16384, 65536, 131072, 524288, 1048576, 36864, 73728, 294912, 589824, 2359296, 122289, 98304, 196608, 786432, 1572864, 13056, 52224, 104448, 417792, 1671168};
//    int p[20] = {2,2,2,2,2,6,6,6,6,6,10,10,10,10,10,20,20,20,20,20};
//    int s[20] = {16,16,16,16,16,72,72,72,72,72,192,192,192,192,192,816,816,816,816,816};
//    N_calc = 1000;
//
//    for(int i=0; i<40; ++i){
//        int ii=i;
//        if (ii>19) ii -= 20;
//        N = n_ptc[ii];
//        n_ptc_box = s[ii];
//        n_rank = p[ii];
//
//        if(i<20){
////////            	 Generate random particle positions.
//            srand (time(NULL));
//            for(unsigned long int i=0;i<N;++i){
//                x[i] = rand();
//                y[i] = rand();
//                z[i] = rand();
//                q[i] = 1;
//            }
//        }
//        else{
//            ////         normal distribution.
//            ////     obtain a seed from the timer
//            std::default_random_engine generator;
//            generator.seed(time(NULL));
//            std::normal_distribution<double> distribution(0.0,1.0);
//            for(unsigned long int i=0;i<N;++i){
//                x[i] = distribution(generator);
//                y[i] = distribution(generator);
//                z[i] = distribution(generator);
//                q[i] = 1;
//            }
//        }
//        scale(x, N);
//        scale(y, N);
//        scale(z, N);
//
//        std::chrono::steady_clock::time_point start, end;
//
//        start = std::chrono::steady_clock::now();
//        Coulomb(x,y,z,q,N,N_calc,phi_check);
//        end = std::chrono::steady_clock::now();
//        auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//
//        start = std::chrono::steady_clock::now();
//        fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
//        end = std::chrono::steady_clock::now();
//        auto t_phi = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//
//        start = std::chrono::steady_clock::now();
//        Coulomb(x,y,z,q,N,N_calc,Ex_check,Ey_check,Ez_check);
//        end = std::chrono::steady_clock::now();
//        auto t2 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//
//        start = std::chrono::steady_clock::now();
//        fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//        end = std::chrono::steady_clock::now();
//        auto t_E = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//
//        double Err_phi = error(phi,phi_check,N_calc);
//        double Err_E = error(Ex,Ey,Ez,Ex_check,Ey_check,Ez_check,N_calc);
//
//        if(i<20) outfile<<"uniform ";
//        else outfile<<"Guassian ";
//
//        outfile<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<Err_phi<<' '<<t_phi<<' '<<t1*N/N_calc<<' '<<Err_E<<' '<<t_E<<' '<<t2*N/N_calc<<endl;
//
//        cout<<i<<endl;
//
//
//
//    }
//
//    outfile.close();
//



//
//
//
//////===========================
////// Demonstrate O(N) efficiency
//////===========================
//
//
//
//    char filename[50] = "0602_Gaussian_field_find_s_rank_6_O(N).txt";
//    std::ofstream outfile;
//    outfile.open(filename);
//    outfile<<"n_ptc rank time"<<endl;
//
//    char filename2[60] = "0602_Gaussian_field_find_s_rank_6_O(N)_average.txt";
//    std::ofstream outfile2;
//    outfile2.open(filename2);
//    outfile2<<"n_ptc rank time"<<endl;
//
//    double n_ptc[10] = {1e5, 2e5, 4e5, 6e5, 8e5, 1e6, 2e6, 4e6, 8e6, 1e7};
//    n_rank = 6;
//    int s[10] = {128, 64, 128, 128, 128, 64, 64, 128, 64, 64};
//
//
//    int n_repeat = 10;
//
//
//    for(int i=0; i<10; ++i){
//            N = n_ptc[i];
//            n_ptc_box = s[i];
//            double t_record = 0;
//            for(int j=0; j<n_repeat; ++j){
//                //////         normal distribution.
//                //////     obtain a seed from the timer
//                std::default_random_engine generator;
//                generator.seed(time(NULL));
//                std::normal_distribution<double> distribution(0.0,1.0);
//                for(unsigned long int i=0;i<N;++i){
//                    x[i] = distribution(generator);
//                    y[i] = distribution(generator);
//                    z[i] = distribution(generator);
//                    q[i] = 1;
//                }
//
//
//                scale(x, N);
//                scale(y, N);
//                scale(z, N);
//                cout<<N<<" particles initialized!"<<endl;
//
//                auto start = std::chrono::steady_clock::now();
//                fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//                auto end = std::chrono::steady_clock::now();
//                auto t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//                cout<<j<<' '<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//                outfile<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//
//                t_record += t3;
//            }
//
//            cout<<i<<' '<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t_record/n_repeat<<endl;
//            outfile2<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t_record/n_repeat<<endl;
//
//    }
//
//
//    outfile.close();
//    outfile2.close();
//
//


//
//////===========================
////// Find the best s for potential calculation with rand 6
//////===========================
//
//
//
//    char filename[50] = "0601_Gaussian_field_find_s_rank_6.txt";
//    std::ofstream outfile;
//    outfile.open(filename);
//    outfile<<"n_ptc rank time"<<endl;
//
//    double n_ptc[10] = {1e5, 2e5, 4e5, 6e5, 8e5, 1e6, 2e6, 4e6, 8e6, 1e7};
//    n_rank = 6;
//    n_ptc_box = 32;
//
//    double t_rec = 1e6;
//
//    for(int ii=0; ii<10; ++ii){
//        N = n_ptc[ii];
//
////        //	 Generate random particle positions.
////        srand (time(NULL));
////        for(unsigned long int i=0;i<N;++i){
////            x[i] = rand();
////            y[i] = rand();
////            z[i] = rand();
////            q[i] = 1;
////        }
//
////////         normal distribution.
//        //////     obtain a seed from the timer
//        std::default_random_engine generator;
//        generator.seed(time(NULL));
//        std::normal_distribution<double> distribution(0.0,1.0);
//        for(unsigned long int i=0;i<N;++i){
//            x[i] = distribution(generator);
//            y[i] = distribution(generator);
//            z[i] = distribution(generator);
//            q[i] = 1;
//        }
//
//
//        scale(x, N);
//        scale(y, N);
//        scale(z, N);
//        cout<<N<<" particles initialized!"<<endl;
//
//
//
//        if (n_ptc_box>32) n_ptc_box /= 2;
//        if (n_ptc_box>128) n_ptc_box = 128;
//
//        auto start = std::chrono::steady_clock::now();
//        fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//        auto end = std::chrono::steady_clock::now();
//        auto t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//        cout<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//        outfile<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//
//        n_ptc_box *= 2;
//
//        while(t3<1.2*t_rec){
//            t_rec = t3;
//
//            start = std::chrono::steady_clock::now();
//            fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//            end = std::chrono::steady_clock::now();
//            t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//            outfile<<N<<' '<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//
//            n_ptc_box *= 2;
//        }
//        t_rec=1e6;
//        n_ptc_box /= 8;
//
//
//
//        if(n_ptc_box>N){
//            t_rec = 1e6;
//            n_ptc_box = 8;
//            break;
//        }
//    }
//
//    outfile.close();



//////=========================================
////// Speed, but doesn't work, need to find good configuration for each n particle
//////=========================================
//
//
//    char filename[50] = "0529_Uniform_n_speed.txt";
//    std::ofstream outfile;
//    outfile.open(filename);
//
//
//    double n_ptc[10] = {1e5, 2e5, 4e5, 6e5, 8e5, 1e6, 2e6, 4e6, 8e6, 1e7};
//    n_rank = 6;
//    n_ptc_box = 128;
//
//
//    for(int ii=0; ii<10; ++ii){
//
//        N = n_ptc[ii];
//        int N_repeat = 10;
//        int tim_E = 0;
//        int tim_phi = 0;
//        int t_E = 0;
//        int t_phi = 0;
//        int err_E = 0;
//        int err_phi = 0;
//
//        for(int i=0;i<N_repeat;++i){
//                ////  normal distribution.
//            ////////     obtain a seed from the timer
//    //        std::default_random_engine generator;
//    //        generator.seed(time(NULL));
//    //        std::normal_distribution<double> distribution(0.0,1.0);
//    //        for(unsigned long int i=0;i<N;++i){
//    //            x[i] = distribution(generator);
//    //            y[i] = distribution(generator);
//    //            z[i] = distribution(generator);
//    //            q[i] = 1;
//    //        }
//        //
//        //	 Generate random particle positions.
//            srand (time(NULL));
//            for(unsigned long int i=0;i<N;++i){
//                x[i] = rand();
//                y[i] = rand();
//                z[i] = rand();
//                q[i] = 1;
//            }
//            scale(x, N);
//            scale(y, N);
//            scale(z, N);
//            cout<<N<<" particles initialized!"<<endl;
//
//            std::chrono::steady_clock::time_point start, end;
//            cout<<"Calculate the first "<<N_calc<<" particles by the pairwise Coulomb formula!"<<endl;
//            start = std::chrono::steady_clock::now();
//            Coulomb(x,y,z,q,N,N_calc,phi_check);
//            end = std::chrono::steady_clock::now();
//            auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<"Coulomb formula: "<<t1<<' '<<N/N_calc*t1<<endl;
//            tim_phi += t1;
//
//
//            cout<<"Calculate the first "<<N_calc<<" particles by the pairwise Coulomb formula!"<<endl;
//            start = std::chrono::steady_clock::now();
//            Coulomb(x,y,z,q,N,N_calc,Ex_check,Ey_check,Ez_check);
//            end = std::chrono::steady_clock::now();
//            t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<"Coulomb formula: "<<t1<<' '<<N/N_calc*t1<<endl;
//            tim_E += t1;
//
//            start = std::chrono::steady_clock::now();
//            fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
//            end = std::chrono::steady_clock::now();
//            auto t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//            cout<<"FMM/Coulomb: "<<t3/(N/N_calc*t1)<<endl;
//
//            double tmp = error(phi_check,phi,N_calc);
//            t_phi += t3;
//            err_phi += tmp;
//
//
//            start = std::chrono::steady_clock::now();
//            fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//            end = std::chrono::steady_clock::now();
//            t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//            cout<<"FMM/Coulomb: "<<t3/(N/N_calc*t1)<<endl;
//            tmp = error(Ex_check,Ey_check,Ez_check,Ex,Ey,Ez,N_calc);
//            outfile<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3<<' '<<tmp<<' '<<t1<<' '<<t3/(N/N_calc*t1)<<endl;
//            t_E += t3;
//            err_E += tmp;
//
//
//        }
//
//        outfile<<n_ptc[ii]<<' '<<t_E/N_repeat<<' '<<t_phi/N_repeat<<' '<<tim_E/N_repeat<<' '<<tim_phi/N_repeat<<' '<<err_E/N_repeat<<' '<<err_phi/N_repeat<<endl;
//    }
//
//    outfile.close();
//
//




//////
//////==================================
////// Calculate field and average error
//////==================================
//
//    char filename[50] = "0528_Potential_Uniform.txt";
//    std::ofstream outfile;
//    outfile.open(filename);
//
//
//
//    int ptc_box[9] = {64,64,64,64,64,64,64,64,64};
//    int order_fmm[9] = {2,3,4,5,6,7,8,9,10};
//    double t[9] = {0};
//    double err[9] = {0};
//    double tim = 0;
//    int N_order = 9;
//    int N_repeat = 10;
//
//    cout<<"start FMM"<<endl;
//
//    for(int i=0;i<N_repeat;++i){
//
//        ////  normal distribution.
//        ////////     obtain a seed from the timer
////        std::default_random_engine generator;
////        generator.seed(time(NULL));
////        std::normal_distribution<double> distribution(0.0,1.0);
////        for(unsigned long int i=0;i<N;++i){
////            x[i] = distribution(generator);
////            y[i] = distribution(generator);
////            z[i] = distribution(generator);
////            q[i] = 1;
////        }
//    //
//    //	 Generate random particle positions.
//        srand (time(NULL));
//        for(unsigned long int i=0;i<N;++i){
//            x[i] = rand();
//            y[i] = rand();
//            z[i] = rand();
//            q[i] = 1;
//        }
//        scale(x, N);
//        scale(y, N);
//        scale(z, N);
//        cout<<N<<" particles initialized!"<<endl;
//
//        std::chrono::steady_clock::time_point start, end;
//        cout<<"Calculate the first "<<N_calc<<" particles by the pairwise Coulomb formula!"<<endl;
//        start = std::chrono::steady_clock::now();
////        Coulomb(x,y,z,q,N,N_calc,Ex_check,Ey_check,Ez_check);
//        Coulomb(x,y,z,q,N,N_calc,phi_check);
//        end = std::chrono::steady_clock::now();
//        auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//        cout<<"Coulomb formula: "<<t1<<' '<<N/N_calc*t1<<endl;
////        outfile<<"Coulomb formula: "<<t1<<' '<<N/N_calc*t1<<endl;
//        tim += t1;
//
//
//        for(int j=0;j<N_order;++j){
//            n_ptc_box = ptc_box[j];
//            n_rank = order_fmm[j];
//            start = std::chrono::steady_clock::now();
////            fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//            fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
//            end = std::chrono::steady_clock::now();
//            auto t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//            cout<<"FMM/Coulomb: "<<t3/(N/N_calc*t1)<<endl;
////            double tmp = error(Ex_check,Ey_check,Ez_check,Ex,Ey,Ez,N_calc);
//            double tmp = error(phi_check,phi,N_calc);
//            outfile<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3<<' '<<tmp<<' '<<t1<<' '<<t3/(N/N_calc*t1)<<endl;
//            t[j] += t3;
//            err[j] += tmp;
//        }
//    }
//    cout<<"end FMM"<<endl;
//    outfile.close();
//
//    outfile.open("0528_Potential_Uniform_avg.txt");
//    for(int i=0; i<N_order; ++i){
//        outfile<<order_fmm[i]<<' '<<ptc_box[i]<<' '<<t[i]/10<<' '<<err[i]/10<<' '<<t[i]/(N/N_calc*tim)<<endl;
//    }
//    outfile.close();



////
//
//    outfile.open("normal_field.txt");
//    std::chrono::steady_clock::time_point start, end;
//
//    cout<<"Calculate the first "<<N_calc<<" particles by the pairwise Coulomb formula!"<<endl;
//    start = std::chrono::steady_clock::now();
//    Coulomb(x,y,z,q,N,N_calc,Ex_check, Ey_check, Ez_check);
//    end = std::chrono::steady_clock::now();
//    auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//    cout<<"Coulomb formula: "<<t1<<' '<<N/N_calc*t1<<endl;
//    outfile<<"Coulomb formula: "<<t1<<' '<<N/N_calc*t1<<endl;
//
//    int ptc_box[6] = {40,60,100,200,400};
//    int order_fmm[4] = {4,6,8};
//
//    cout<<"start FMM"<<endl;
//
//    for(int i=0;i<5;++i){
//        for(int j=0;j<3;++j){
//            n_ptc_box = ptc_box[i];
//            n_rank = order_fmm[j];
//            start = std::chrono::steady_clock::now();
//            fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//            end = std::chrono::steady_clock::now();
//            auto t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//            cout<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3<<endl;
//            cout<<"FMM/Coulomb: "<<t3/(N/N_calc*t1)<<endl;
//            outfile<<"FMM: "<<n_rank<<' '<<n_ptc_box<<' '<<t3<<' '<<error(Ex_check,Ey_check,Ez_check,Ex,Ey,Ez,N_calc)<<endl;
//        }
//    }
//    cout<<"end FMM"<<endl;
//    outfile.close();

//
//    cout<<"Start FMM for field ... "<<endl;
//    fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, Ex, Ey, Ez);
//    cout<<"End FMM!"<<endl;
//    Coulomb(x,y,z,q,N,N_calc,Ex_check, Ey_check, Ez_check);
//    cout<<"Error: "<<error(Ex_check,Ey_check,Ez_check,Ex,Ey,Ez,N_calc)<<endl;
//
//    cout<<"Start FMM for potential... "<<endl;
//    fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
//    cout<<"End FMM!"<<endl;
//    Coulomb(x,y,z,q,N,N_calc,phi_check);
//    cout<<"Error: "<<error(phi_check,phi,N_calc)<<endl;

//    char filename[30] = "check_field.txt";
//    std::ofstream outfile;
//    outfile.open(filename);
//    for(unsigned long int i=0;i<N_calc;++i){
//        outfile<<Ex[i]<<' '<<Ey[i]<<' '<<Ez[i]<<' '<<Ex_check[i]<<' '<<Ey_check[i]<<' '<<Ez_check[i]<<' '<<endl;
//    }
//    outfile.close();

//    return 0;

//    std::chrono::steady_clock::time_point start, end;
//
//    cout<<"Calculate the first "<<N_calc<<" particles by the pairwise Coulomb formula!"<<endl;
//    start = std::chrono::steady_clock::now();
//    Coulomb(x,y,z,q,N,N_calc,phi_check);
//    end = std::chrono::steady_clock::now();
//    auto t1 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//    cout<<"Coulomb formula: "<<t1<<' '<<N/N_calc*t1<<endl;
//
//    cout<<"start FMM"<<endl;
//    start = std::chrono::steady_clock::now();
//    fmm(x,  y,  z,  q, N, n_rank, n_ptc_box, phi);
//    end = std::chrono::steady_clock::now();
//    auto t3 = 0.001*std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//    cout<<"FMM: "<<t3<<endl;
//    cout<<"end FMM"<<endl;
//    cout<<"FMM/Coulomb: "<<t3/(N/N_calc*t1)<<endl;

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

//    cout<<"Relative error: "<<error(phi,phi_check,N_calc)<<endl;

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
