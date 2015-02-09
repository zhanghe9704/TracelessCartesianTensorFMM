/**********************************
Function.cpp
Declare the math functions used in the kernel functions for the multiple level fast multipole algorithm using tensors

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/


#ifndef FUNCTION_HPP
#define FUNCTION_HPP

int Find_index(int n1, int n2, int n3);
double combination(int n, int m);
double combination_HH(int n, int m);
void Symmetric_Tensor(double x, double y, double z, double *SymmeticTensor);
double Nabla_1_element_r(int n1, int n2, int n3, int n, double x, double y, double z, double r_2, double r_coe);
void Nabla_r_direct(double end_x, double end_y, double end_z, double begin_x, double begin_y, double begin_z, double *Nabla_R);
void Nabla_r_traceless(double end_x, double end_y, double end_z, double begin_x, double begin_y, double begin_z, double *Nabla_R);
void Contraction(double *High_rank_Tensor, double *Low_rank_Tensor, double *HL_rank_Tensor, int m, int n);
void Contraction_traceless(double *High_rank_Tensor, double *Low_rank_Tensor, double *HL_rank_Tensor, int m, int n);

#endif
