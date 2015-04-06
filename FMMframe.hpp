/**********************************
FMMframe.hpp

Declare the function to calculate the Coulomb potential using MLFMA

version 2.0
By He Zhang, 04/06/2015
Calculate 3D field

version 1.0
By He Zhang, 02/2015

***********************************/

#ifndef FMMFRAME_HPP
#define FMMFRAME_HPP

//Calculate the Coulomb field on each particle by MLFMA, NOT finished
int fmm(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * Ex, double * Ey, double * Ez);

//Calculate the Coulomb potential on each particle by MLFMA
int fmm(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * phi);

int fmm(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * phi, double * Ex, double * Ey, double * Ez);

int fmm_dr(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * Ex, double * Ey, double * Ez);

#endif
