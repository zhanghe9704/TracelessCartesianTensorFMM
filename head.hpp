/**********************************
head.hpp
List of the head files and global variables

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/


#ifndef HEAD_HPP
#define HEAD_HPP

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "FMMConstants.hpp"
#include "Function.hpp"
#include "FMMKernel.hpp"


#include "FMMframe.hpp"
#include "box.hpp"

extern int n_Max_rank, Number_of_total_element;
extern unsigned long int Number_of_particle;
extern double *scratch;
extern double *scratch2;
extern unsigned long int *ptclist;
extern double * multipole_expns;
extern double * local_expns;


#endif
