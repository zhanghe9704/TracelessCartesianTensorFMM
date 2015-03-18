/**********************************
Function.cpp
Declare the math functions used in the kernel functions for the multiple level fast multipole algorithm using tensors


version 2.0
By He Zhang 03/06/2015
All functions revised and more functions added for better performance.

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/


#ifndef FUNCTION_HPP
#define FUNCTION_HPP

//output i3>=i2>=i1, if sequence changed return a number greater than 0;
int sequence3 (int x[3], int idx[3]);

//Calculate the totally symmetric tensor r^n
void Symmetric_Tensor(double x, double y, double z, double *SymmeticTensor);

//Calculate the value of a specific element of the tensor/operator Nabla 1/r
double Nabla_1_element_r(int n1, int n2, int n3, int n, double x, double y, double z, double r_2, double r_coe, int &cnt, double * coef);

//Calculate the tensor Nabla 1/r for a give r(x,y,z)
void Nabla_r_traceless(double x, double y, double z, double * coef, double *Nabla_R);

//Contraction of two tensors. The high rank one is traceless totally symmetric, the low rank one is totally symmetric, the result is traceless totally symmetric
void Contraction_traceless(double *High_rank_Tensor, double *Low_rank_Tensor, double *HL_rank_Tensor, int m, int n);

//Using the value of the first 2*n+1 elements to calculate the value of the other elements for a traceless totally symmetric tensor
void fill_traceless_tensor(double * tensor);

//fill the totally symmetric tensor r^n
void fill_symmetric_tensor_r(double r2, double *SymmetricTensor);
//contraction of two tensors of the same rank
double Contraction_equal_rank(double * Tensor1, double * Tensor2, int n);

//Detracer operator, convert a totally symmetric tensor into a traceless totally symmetric tensor. Not used in current version.
void detracer(double * symmetric_tensor, double * traceless_tensor);


//Take derivative before contraction of two tensors of the same rank
void Contraction_dr(double * Tensor1, double * Tensor2, int n, double x, double y, double z, double &cx, double &cy, double &cz);

//Take derivative of the tensor/operator Nabla 1/r for a given r
void Nabla_1_element_r_dr(int n1, int n2, int n3, int n, double x, double y, double z, double r_2, double r_coe, int &cnt, double * coef, double * one_element);

//Calculate the derivative of the tensor Nabla 1/r for a give r(x,y,z)
void Nabla_r_dr(double x, double y, double z, double * coef, double *Nabla_x, double *Nabla_y, double *Nabla_z);
#endif
