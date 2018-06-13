#ifndef GLOBAL_H
#define GLOBAL_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>

#include "box.h"
#include "fmm_constants.h"
#include "fmm_frame.h"
#include "fmm_kernel.h"
#include "functions.h"

extern int MAX_RANK;
extern int TOTAL_ELEMENT_NUMBER;
extern unsigned long int TOTAL_PARTICLE_NUMBER;
extern double *g_scratch;
extern double *g_scratch2;
extern unsigned long int *g_ptclist;
extern double *g_multipole_expns;
extern double *g_local_expns;
extern double *g_combination_coef;
extern double *g_nabla_coef;
extern double *g_pow_x;
extern double *g_pow_y;
extern double *g_pow_z;
extern double *g_pow_r2;
extern Flag g_flag;

#endif
