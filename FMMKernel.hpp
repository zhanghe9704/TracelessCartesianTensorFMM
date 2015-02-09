/**********************************
FMMKernel.hpp
Declare the Kernel functions for the multiple level fast multipole algorithm using tensors

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/



#ifndef FMMKERNEL_HPP
#define FMMKERNEL_HPP

#include "box.hpp"

// void Charge_to_Multipole(int Charge_Number, double *Charge, double Multipole_x, double Multipole_y, double Multipole_z, double *Charge_x, double *Charge_y, double *Charge_z, double *Born_Multipole, Config& config);
void Charge_to_Multipole(int Charge_Number, double *Charge, double Multipole_x, double Multipole_y, double Multipole_z, double *Charge_x, double *Charge_y, double *Charge_z, double *Born_Multipole);
void Charge_to_Multipole(Box & box, double *Charge, double *Charge_x, double *Charge_y, double *Charge_z, unsigned long int *ptclist, double *Born_Multipole);


void Multipole_to_Multipole(double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double *Old_M, double *New_M);

void Multipole_to_Local(double *Multipole_for_trans, double Multi_x, double Multi_y, double Multi_z, double Local_x, double Local_y, double Local_z, double *M2L_translation);

void Local_to_Local(double *Old_Local, double Old_x, double Old_y, double Old_z, double New_x, double New_y, double New_z, double *New_Local);

double MultipolePotential(double *Multipole, double Multi_x, double Multi_y, double Multi_z, double Poten_x, double Poten_y, double Poten_z);

void Charge_to_Local(double charge_number, double *q, double *old_x, double *old_y, double *old_z, double new_x, double new_y, double new_z, double *L_expansion);

double LocalPotential(double *Local_expan, double Local_x, double Local_y, double Local_z, double observer_x, double observer_y, double observer_z);

void Charge_to_Multipole(Box & box, double *Charge, double *Charge_x, double *Charge_y, double *Charge_z, double *ptclist, double *Born_Multipole);
#endif
