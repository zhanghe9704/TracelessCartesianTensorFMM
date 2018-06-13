#ifndef FMM_FRAME_H
#define FMM_FRAME_H

#include "fmm_frame.h"

//Calculate both the potential and the field by MLFMA
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, unsigned int max_rank, unsigned int n_ptc_box,
        double *phi, double *ex, double *ey, double *ez);

//Calculate the Coulomb field on each particle by MLFMA
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, unsigned int max_rank, unsigned int n_ptc_box,
        double *ex, double *ey, double *ez);

//Calculate the Coulomb potential on each particle by MLFMA
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, unsigned int max_rank, unsigned int n_ptc_box,
        double *phi);
#endif
