/**********************************
FMMframe.hpp

Declare the function to calculate the Coulomb potential using MLFMA

version 1.0
By He Zhang, 02/2015

***********************************/

#ifndef FMMFRAME_HPP
#define FMMFRAME_HPP

//Calculate the Coulomb field on each particle by MLFMA, NOT finished
int fmm(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * Ex, double * Ey, double * Ez);

//Calculate the Coulomb potential on each particle by MLFMA
int fmm(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * phi);

#endif
