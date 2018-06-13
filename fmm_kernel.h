#ifndef FMM_KERNEL_H
#define FMM_KERNEL_H

#include "box.h"

//Calculate the multipole from charges inside a childless box
void charge_to_multipole(Box &box, double *charge, double *charge_x, double *charge_y, double *charge_z,
                         double *multipole);
//Translate multipole from old position to new position
void multipole_to_multipole(double old_x, double old_y, double old_z, double new_x, double new_y, double new_z,
                            double *multipole_coef, double *old_mulitpole, double *new_multipole);
//Calculate the coef for the M2M translation operator
void multipole_to_multipole_coef(double boxsize, double * multipole_coef);
//Update the M2M operator for one level up (go upwards)
void update_multipole_to_multipole_coef(double * multipole_coef);

//Calculate the M2L operator nabla_r, r is the vector that links the centers of the two boxes
void calc_nabla_r(double boxsize, double *nabla_r);
//Update the M2L operator for one level down (go downwards)
void update_nabla_r(double *nabla_r);
//Convert a multipole expansion into a local expansion
void multipole_to_local(double *multipole, double multi_x, double multi_y, double multi_z, double local_x,
                        double local_y, double local_z, double *saved_nabla_r, double boxsize, double *local_expn);
//Translate the local expansion from a parent box to a child box
void local_to_local(double *old_local, double old_x, double old_y, double old_z, double new_x, double new_y,
                    double new_z, double *origin_rho_tensor, double *new_local);
//Calc the L2L operator
void calc_rho_tensor(double boxsize, double *rho_tensor);
//Update the L2L operator for one level down (go downwards)
void update_rho_tensor(double *rho_tensor);
//Calculate the local expansion from charges inside a ill-separated box
void charge_to_local_traceless(Box &box, double *q, double *old_x, double *old_y, double *old_z, double new_x,
                               double new_y, double new_z, double *local_expn);
//Evaluate the kernel function on charges by the multipole expansion
int multipole_to_charge(double *multipole, double multi_x, double multi_y, double multi_z, double obsv_x, double obsv_y,
                        double obsv_z, double &phi, double &ex, double &ey, double &ez);
//Evaluate the kernel function on charges by the local expansion
int local_to_charge(double *local_expn, double local_x, double local_y, double local_z, double obsv_x, double obsv_y,
                    double obsv_z, double &phi, double &ex, double &ey, double &ez);
#endif
