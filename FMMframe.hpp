#ifndef FMMFRAME_HPP
#define FMMFRAME_HPP

int fmm(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * Ex, double * Ey, double * Ez);
int fmm(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * phi);

#endif
