#include "fmm_constants.h"

//Global variables
int MAX_RANK;
int TOTAL_ELEMENT_NUMBER;
unsigned long int TOTAL_PARTICLE_NUMBER;
double * g_scratch;
double * g_scratch2;
unsigned long int * g_ptclist;
double * g_multipole_expns;
double * g_local_expns;
double * g_combination_coef;
double * g_nabla_coef;
double * g_pow_x;
double * g_pow_y;
double * g_pow_z;
double * g_pow_r2;
//flag indicates what's to calculate.
Flag g_flag = Flag::BOTH;

//Calculate the coefs for kCombination operator
int calc_combination_coef(double * g_combination_coef) {
    for(int i=0; i<TOTAL_ELEMENT_NUMBER; ++i){
        int n1, n2, n3;
        n1 = kIndexN1[i];
        n2 = kIndexN2[i];
        n3 = kIndexN3[i];
        g_combination_coef[i] = kFactorial[n1+n2+n3] * kInvFactorial[n1] * kInvFactorial[n2] * kInvFactorial[n3];
    }
    return 0;
}

//Calculate the coefs for Nabla operator
int calc_nabla_1_emement_coef(double * g_nabla_coef) {
    int cnt = 0;
    for(int n=0; n<MAX_RANK+1;++n) {
		for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i) {
			int n1, n2, n3;
			n1 = kIndexN1[i];
            n2 = kIndexN2[i];
            n3 = kIndexN3[i];
            for (int m1 = 0; m1 <= (n1 / 2); m1++) {
                for (int m2 = 0; m2 <= (n2 / 2); m2++) {
                    for (int m3 = 0; m3 <= (n3 / 2); m3++) {
                        int m = m1 + m2 + m3;
                        g_nabla_coef[cnt] = kPowerNegOne[m] * kCombinationHH[n1][m1] * kCombinationHH[n2][m2] *
                                            kCombinationHH[n3][m3] * kFactorialOdd[n - m] ;
                        ++cnt;
                    }
                }
            }
		}
	}
	return 0;
}

//Set the value for global invariables and initialize g_scratching arrays
int configure_fmm(int max_rank, unsigned long int n_ptc, unsigned long int n_box){

	MAX_RANK = max_rank;
	TOTAL_ELEMENT_NUMBER = kRankNTensorStart[MAX_RANK+1];
	TOTAL_PARTICLE_NUMBER = n_ptc;
	//Create two g_scratch arrays in heap
	g_scratch = new double[TOTAL_ELEMENT_NUMBER];
	memset(g_scratch, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));
	g_scratch2 = new double[TOTAL_ELEMENT_NUMBER];
	memset(g_scratch2, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));

    // Create an array to store all the multipole expansions of all the boxes
	// the starting address for the i-th multipole is &multipole[i*TOTAL_ELEMENT_NUMBER], i counts from zero.
	g_multipole_expns = new double[n_box*TOTAL_ELEMENT_NUMBER];
	memset(g_multipole_expns, 0, n_box*TOTAL_ELEMENT_NUMBER*sizeof(double));

    // Create an array to store all the local expansions of all the boxes
	// the starting address for the i-th local is &multipole[i*TOTAL_ELEMENT_NUMBER], i counts from zero.
	g_local_expns = new double[n_box*TOTAL_ELEMENT_NUMBER];
	memset(g_local_expns, 0, n_box*TOTAL_ELEMENT_NUMBER*sizeof(double));

    //Save the coefs for kCombination operator
    g_combination_coef = new double[TOTAL_ELEMENT_NUMBER];
    calc_combination_coef(g_combination_coef);

    //Use to save the power of x, y, z, r2 up to the max rank
    g_pow_x = new double[MAX_RANK+1];
    g_pow_y = new double[MAX_RANK+1];
    g_pow_z = new double[MAX_RANK+1];
    g_pow_r2 = new double[MAX_RANK+1];
    memset(g_pow_x, 0, (MAX_RANK+1)*sizeof(double));
    memset(g_pow_y, 0, (MAX_RANK+1)*sizeof(double));
    memset(g_pow_z, 0, (MAX_RANK+1)*sizeof(double));
    memset(g_pow_r2, 0, (MAX_RANK+1)*sizeof(double));

	return 0;
}

//Release the g_scratch memory in heap
int end_fmm(){
	delete[] g_scratch;
	delete[] g_multipole_expns;
	delete[] g_local_expns;
	delete[] g_scratch2;
	delete[] g_combination_coef;
	delete[] g_pow_x;
	delete[] g_pow_y;
	delete[] g_pow_z;
	delete[] g_pow_r2;

	return 0;
}
