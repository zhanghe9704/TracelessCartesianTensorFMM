/**********************************
FMMConstants.cpp
Define and initialize the constants used in the kernel functions for the multiple level fast multipole algorithm using tensors

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/

#include "FMMConstants.hpp"

//Global variables
int n_Max_rank, Number_of_total_element;
unsigned long int Number_of_particle;
double * scratch;
double * scratch2;
unsigned long int * ptclist;
double * multipole_expns;
double * local_expns;
double * combination_coef;
double * Nabla_1_element_r_coef;
double * pow_x;
double * pow_y;
double * pow_z;
double * pow_r2;

//int nabla_idx[4][4][2]={0};

int calc_combination_coef(double * combination_coef){
    for(int i=0; i<Number_of_total_element; ++i){
        int n1, n2, n3;
        n1 = index_n1[i];
        n2 = index_n2[i];
        n3 = index_n3[i];

        combination_coef[i] = Factorial[n1+n2+n3] * inv_Factorial[n1] * inv_Factorial[n2] * inv_Factorial[n3];
    }
    return 0;
}

int Calc_Nabla_1_emement_coef(double * Nabla_1_element_r_coef){
    int cnt = 0;
    for(int n=0; n<n_Max_rank+1;++n){

		for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
			int n1, n2, n3;
			n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            for (int m1 = 0; m1 <= (n1 / 2); m1++)
                {
                for (int m2 = 0; m2 <= (n2 / 2); m2++)
                {
                    for (int m3 = 0; m3 <= (n3 / 2); m3++)
                    {
                        int m = m1 + m2 + m3;
                        Nabla_1_element_r_coef[cnt] = order_minus_one[m] * combination_HH[n1][m1] * combination_HH[n2][m2] *combination_HH[n3][m3] * Factorial_odd[n - m] ;
                        ++cnt;
                    }
                }
            }
		}
	}
	return 0;
}
/*
 * int calc_nabla_idx(){
 *     int shift = 0;
 *     for(int i3=2; i3<4; ++i3){
 *         for(int i2=0; i2<i3+1; ++i2){
 *             for(int i1=0; i1<i2+1; ++i1){
 *                 nabla_idx[i1][i2][i3-2] = shift;
 *                 ++shift;
 *             }
 *         }
 *     }
 *     return 0;
 * }
 */

int configure_fmm(int Max_rank, unsigned long int n_ptc, unsigned long int n_box){

	n_Max_rank = Max_rank;
	Number_of_total_element = n_Rank_Multipole_Start_Position[n_Max_rank+1];
	Number_of_particle = n_ptc;
	//Create two scratch arrays in heap
	scratch = new double[Number_of_total_element];
	memset(scratch, 0, Number_of_total_element*sizeof(double));
	scratch2 = new double[Number_of_total_element];
	memset(scratch2, 0, Number_of_total_element*sizeof(double));

    // Create an array to store all the multipole expansions of all the boxes
	// the starting address for the i-th multipole is &multipole[i*Number_of_total_element], i counts from zero.
	multipole_expns = new double[n_box*Number_of_total_element];
	memset(multipole_expns, 0, n_box*Number_of_total_element*sizeof(double));

	local_expns = new double[n_box*Number_of_total_element];
	memset(local_expns, 0, n_box*Number_of_total_element*sizeof(double));

    combination_coef = new double[Number_of_total_element];
    calc_combination_coef(combination_coef);

    Nabla_1_element_r_coef = new double[Nabla_1_element_r_length[Max_rank]];
    Calc_Nabla_1_emement_coef(Nabla_1_element_r_coef);

    pow_x = new double[Max_rank+1];
    pow_y = new double[Max_rank+1];
    pow_z = new double[Max_rank+1];
    pow_r2 = new double[Max_rank+1];
    memset(pow_x, 0, (Max_rank+1)*sizeof(double));
    memset(pow_y, 0, (Max_rank+1)*sizeof(double));
    memset(pow_z, 0, (Max_rank+1)*sizeof(double));
    memset(pow_r2, 0, (Max_rank+1)*sizeof(double));

	return 0;
}

//Release the scratch memory in heap
int end_fmm(){
	delete[] scratch;
	delete[] multipole_expns;
	delete[] local_expns;
	delete[] scratch2;
	delete[] combination_coef;
	delete[] Nabla_1_element_r_coef;
	delete[] pow_x;
	delete[] pow_y;
	delete[] pow_z;
	delete[] pow_r2;

	return 0;
}
