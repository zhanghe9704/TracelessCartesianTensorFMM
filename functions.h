#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//output i3>=i2>=i1, if sequence changed return a number greater than 0;
int sequence3 (int x[3], int idx[3]);

//Calculate the totally symmetric tensor r^n
void symmetric_tensor(double x, double y, double z, double *tensor);

//Calculate the totally symmetric tensor r^n times charge Q for charge to multipole calculation
//void Symmetric_Tensor_C2M(double q, double x, double y, double z, double *SymmetricTensor);
void symmetric_tensor_c2m(double q, double x, double y, double z, double *symmetric_tensor);

//Calculate the value of a specific element of the tensor/operator Nabla 1/r
double nabla_element_r(int n1, int n2, int n3, int n, double x, double y, double z, double r_2, double r_coe, int &cnt,
                       double *coef);

//Calculate the tensor Nabla 1/r for a give r(x,y,z)
void nabla_r_traceless(double x, double y, double z, double *coef, double *nabla_r);

//Using the value of the first 2*n+1 elements to calculate the value of the other elements for a traceless totally
//symmetric tensor
void fill_traceless_tensor(double *tensor);

//fill the totally symmetric tensor r^n
void fill_symmetric_tensor_r(double r2, double *symmetric_tensor);
//contraction of two tensors of the same rank
double contraction_equal_rank(double * tensor1, double * tensor2, int n);

//Contraction of two tensors. The high rank one is traceless totally symmetric, the low rank one is totally symmetric,
//the result is traceless totally symmetric
void contraction_traceless(double *high_rank_tensor, double *low_rank_tensor, double *output_tensor, int m, int n);

//Take derivative before contraction of two tensors of the same rank
void contraction_dr(double *tensor1, double *tensor2, int n, double x, double y, double z, double &cx, double &cy,
                    double &cz);

#endif
