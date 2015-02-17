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

int configure_fmm(int Max_rank, unsigned long int n_ptc, unsigned long int n_box){

	n_Max_rank = Max_rank;
	Number_of_total_element = n_Rank_Multipole_Start_Position[n_Max_rank+1];
	Number_of_particle = n_ptc;
	//Create two scratch arrays in heap
	scratch = new double[Number_of_total_element];
	memset(scratch, 0, Number_of_total_element*sizeof(double));
	scratch2 = new double[Number_of_total_element];
	memset(scratch2, 0, Number_of_total_element*sizeof(double));

//	ptclist = new unsigned long int[n_ptc];
//	memset(ptclist, 0, n_ptc*sizeof(unsigned long int));

    // Create an array to store all the multipole expansions of all the boxes
	// the starting address for the i-th multipole is &multipole[i*Number_of_total_element], i counts from zero.
	multipole_expns = new double[n_box*Number_of_total_element];
	memset(multipole_expns, 0, n_box*Number_of_total_element*sizeof(double));

	local_expns = new double[n_box*Number_of_total_element];
	memset(local_expns, 0, n_box*Number_of_total_element*sizeof(double));


	return 0;
}

//Release the scratch memory in heap
int end_fmm(){
	delete[] scratch;
	delete[] multipole_expns;
	delete[] local_expns;
	delete[] scratch2;
	return 0;
}
