/**********************************
FMMKernel.hpp
Declare the Kernel functions for the multiple level fast multipole algorithm using tensors

version 3.0
By He Zhang & He Huang, 04/06/2015
Calculate 3D field.

version 2.0
By He Zhang, 03/06/2015
All functions revised for better performance.

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/



#ifndef FMMKERNEL_HPP
#define FMMKERNEL_HPP

#include "box.hpp"
#include "head.hpp"

//Calculate the multipole from charges inside a childless box
void Charge_to_Multipole(Box & box, double *Charge, double *Charge_x, double *Charge_y, double *Charge_z, unsigned long int *ptclist, double *Born_Multipole);

//Translate multipole from old position to new position
void Multipole_to_Multipole(double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double * multipole_coef, double *Old_M, double *New_M);
//Calculate the coef for the M2M translation operator
void multipole_to_multipole_coef(double boxsize, double * multipole_coef);
//Update the M2M operator for one level up (go upwards)
void update_multipole_to_multipole_coef(double * multipole_coef);

//Calculate the M2L operator Nabla_R, R is the vector that links the centers of the two boxes
void Calc_Nabla_R(double boxsize, double * Nabla_R);
//Update the M2L operator for one level down (go downwards)
void update_Nabla_R(double * Nabla_R);
//Convert a multipole expansion into a local expansion
void Multipole_to_Local(double *Multipole_for_trans, double Multi_x, double Multi_y, double Multi_z, double Local_x, double Local_y, double Local_z, double * Saved_Nabla_R, double boxsize, double *M2L_translation);

//Translate the local expansion from a parent box to a child box
void Local_to_Local(double *Old_Local, double Old_x, double Old_y, double Old_z, double New_x, double New_y, double New_z, double * lOrigin_Rho_Tensor, double *New_Local);
//Calc the L2L operator
void Calc_Rho_Tensor(double boxsize, double * Rho_Tensor);
//Update the L2L operator for one level down (go downwards)
void update_Rho_Tensor(double * Rho_Tensor);

//Calculate the Coulomb potential using the multipole expansion, for ill-separated boxes
double MultipolePotential(double *Multipole, double Multi_x, double Multi_y, double Multi_z, double Poten_x, double Poten_y, double Poten_z);

//Calculate the local expansion from charges inside a ill-separated box
void Charge_to_Local_traceless(Box &box, unsigned long int * ptclist, double *q, double *old_x, double *old_y, double *old_z, double new_x, double new_y, double new_z, double *L_expansion);

//Calculate the Coulomb potential from the local expansion of a childless box
double LocalPotential(double *Local_expan, double Local_x, double Local_y, double Local_z, double observer_x, double observer_y, double observer_z);

void MultipoleField(double *Multipole, double Multi_x, double Multi_y, double Multi_z, double Poten_x, double Poten_y, double Poten_z, double &Ex, double &Ey, double &Ez);

void LocalField(double *Local_expan, double Local_x, double Local_y, double Local_z, double observer_x, double observer_y, double observer_z, double &Ex, double &Ey, double &Ez);

#endif
