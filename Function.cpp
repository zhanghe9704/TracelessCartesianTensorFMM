/**********************************
Function.cpp
Define the math functions used in the kernel functions for the multiple level fast multipole algorithm using tensors

version 3.0
By He Zhang 04/06/2015
Add functions needed for field calculation.

version 2.0
By He Zhang 03/06/2015
All functions revised and more functions added for better performance.

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/


#include "head.hpp"

//output i3>=i2>=i1, if sequence changed return a number greater than 0;
int sequence3 (int x[3], int idx[3]){
    int change = 0;
    for (int i=0; i<3; ++i) idx[i] = i;
    if (x[1]>x[2])  {idx[2]=1; idx[1]=2; ++change; }
    if (x[0]>x[idx[2]]) {idx[0]=idx[2]; idx[2]=0; ++change;}
    if (x[idx[0]]>x[idx[1]]) { int tmp = idx[0]; idx[0] = idx[1]; idx[1] = tmp; ++change;}
    return change;
}

//fill the totally symmetric tensor r^n
void fill_symmetric_tensor_r(double r2, double *SymmetricTensor){
    for(int n=0; n<n_Max_rank+1; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]+2*n+1; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
            int n1, n2, n3, index_a, index_b, index;
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            index_a = find_index[n1 + 2][n2][n3 - 2];
            index_b = find_index[n1][n2 + 2][n3 - 2];
            index = find_index[n1][n2][n3-2];
            SymmetricTensor[i] = SymmetricTensor[index]*r2-SymmetricTensor[index_a]-SymmetricTensor[index_b];
        }
    }
}

//Calculate the totally symmetric tensor r^n
void Symmetric_Tensor(double x, double y, double z, double *SymmetricTensor){
    SymmetricTensor[0] = 1;
    if (n_Max_rank>0){
        SymmetricTensor[1] = x;
        SymmetricTensor[2] = y;
        SymmetricTensor[3] = z;
    }
    double r2 = x*x+y*y+z*z;
    for(int n=2; n<n_Max_rank+1; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
            int n1 = index_n1[i];
            int n2 = index_n2[i];
            int n3 = index_n3[i];
            if (n1>0){
                int index = find_index[n1-1][n2][n3];
                SymmetricTensor[i] = SymmetricTensor[index]*x;
            }
            else{
                int index = find_index[n1][n2-1][n3];
                SymmetricTensor[i] = SymmetricTensor[index]*y;
            }
        }
    }

    fill_symmetric_tensor_r(r2, SymmetricTensor);

}


//Calculate the totally symmetric tensor r^n times charge Q for charge to multipole calculation
void Symmetric_Tensor_C2M(double q, double x, double y, double z, double *SymmetricTensor){
    SymmetricTensor[0] = 1;
    if (n_Max_rank>0){
        SymmetricTensor[1] = x*q;
        SymmetricTensor[2] = y*q;
        SymmetricTensor[3] = z*q;
    }
    double r2 = x*x+y*y+z*z;
    for(int n=2; n<n_Max_rank+1; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
            int n1 = index_n1[i];
            int n2 = index_n2[i];
            int n3 = index_n3[i];
            if (n1>0){
                int index = find_index[n1-1][n2][n3];
                SymmetricTensor[i] = SymmetricTensor[index]*x;
            }
            else{
                int index = find_index[n1][n2-1][n3];
                SymmetricTensor[i] = SymmetricTensor[index]*y;
            }
        }
    }

    fill_symmetric_tensor_r(r2, SymmetricTensor);

}

//Using the value of the first 2*n+1 elements to calculate the value of the other elements for a traceless totally symmetric tensor
void fill_traceless_tensor(double * tensor){
    for(int n=0; n<=n_Max_rank; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]+2*n+1; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
            int n1, n2, n3, index_a, index_b;
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            index_a = find_index[n1 + 2][n2][n3 - 2];
            index_b = find_index[n1][n2 + 2][n3 - 2];
            tensor[i] = -tensor[index_a] - tensor[index_b];
        }
	}

}

//Calculate the value of a specific element of the tensor/operator Nabla 1/r
double Nabla_1_element_r(int n1, int n2, int n3, int n, double x, double y, double z, double r_2, double r_coe, int &cnt, double * coef)
{
	// reason of inputting "n, r^2, r_coe" is to accelerate calculation for the other fuctions
	// r_2 = x^2 + y^2 + z^2
	// r_coe = (-1)^n * r^(-2n)
	// n = n1+n2+n3

	double one_element = 0;

	for (int m1 = 0; m1 <= (n1 / 2); m1++){
		for (int m2 = 0; m2 <= (n2 / 2); m2++){
			for (int m3 = 0; m3 <= (n3 / 2); m3++){
				int m = m1 + m2 + m3;
				one_element += coef[cnt] * pow_r2[m]*pow_x[n1 - 2 * m1]*pow_y[n2 - 2 * m2]*pow_z[n3 - 2 * m3];
				++cnt;
			}
		}
	}
	one_element *= r_coe;

	return one_element;
}

//Take derivative of the tensor/operator Nabla 1/r for a given r
void Nabla_1_element_r_dr(int n1, int n2, int n3, int n, double x, double y, double z, double r_2, double r_coe, int &cnt, double * coef, double * one_element)
{
	int nn[3] = {n1,n2,n3};
	double c1[3] = {1/x, 1/y, 1/z};
	double c2[3] = {x/r_2, y/r_2, z/r_2};
	int c3 = -2*n-1;

	for(int i=0; i<3; ++i) one_element[i] = 0;

	for (int m1 = 0; m1 <= (n1 / 2); m1++){
		for (int m2 = 0; m2 <= (n2 / 2); m2++){
			for (int m3 = 0; m3 <= (n3 / 2); m3++){
				int m = m1 + m2 + m3;
                int mm[3] = {m1,m2,m3};
				double tmp = coef[cnt] * pow_r2[m]*pow_x[n1 - 2 * m1]*pow_y[n2 - 2 * m2]*pow_z[n3 - 2 * m3];
				for(int id=0; id<3; ++id){
                    double tmp2 = tmp*((nn[id] - 2*mm[id])*c1[id]+(c3+2*m)*c2[id]);
                    one_element[id] += tmp2;
				}
				++cnt;
			}
		}
	}

	for(int i=0; i<3; ++i) one_element[i] *= r_coe;
}


//Calculate the derivative of the tensor Nabla 1/r for a give r(x,y,z)
void Nabla_r_dr(double x, double y, double z, double * coef, double *Nabla_x, double *Nabla_y, double *Nabla_z)
{
    double r_2 = x*x + y*y + z*z;
    double inv_r2 = 1.0/r_2;
    double inv_r = sqrt(inv_r2);

    pow_x[0] = 1;
    pow_y[0] = 1;
    pow_z[0] = 1;
    pow_r2[0] = 1;

    for(int i=1; i<n_Max_rank+1; ++i){
        pow_x[i] = x*pow_x[i-1];
        pow_y[i] = y*pow_y[i-1];
        pow_z[i] = z*pow_z[i-1];
        pow_r2[i] = r_2*pow_r2[i-1];
    }

    //Calculate the first 2n+1 elements for each rank n
    int cnt = 0;
    double one_element[3] = {0,0,0};
	for(int n=0; n<n_Max_rank+1; ++n){
        double r_coe = order_minus_one[n] * pow(inv_r2, n) * inv_r;
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
            int n1, n2, n3;
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            Nabla_1_element_r_dr(n1, n2, n3, n, x, y, z, r_2, r_coe, cnt, coef, one_element);
            Nabla_x[i] = one_element[0];
            Nabla_y[i] = one_element[1];
            Nabla_z[i] = one_element[2];
        }
	}
}



//Calculate the tensor Nabla 1/r for a give r(x,y,z)
void Nabla_r_traceless_formula(double x, double y, double z, double * coef, double *Nabla_R)
{
    double r_2 = x*x + y*y + z*z;
    double inv_r2 = 1.0/r_2;
    double inv_r = sqrt(inv_r2);

    pow_x[0] = 1;
    pow_y[0] = 1;
    pow_z[0] = 1;
    pow_r2[0] = 1;

    for(int i=1; i<n_Max_rank+1; ++i){
        pow_x[i] = x*pow_x[i-1];
        pow_y[i] = y*pow_y[i-1];
        pow_z[i] = z*pow_z[i-1];
        pow_r2[i] = r_2*pow_r2[i-1];
    }

    //Calculate the first 2n+1 elements for each rank n
    int cnt = 0;
	for(int n=0; n<n_Max_rank+1; ++n){
        double r_coe = order_minus_one[n] * pow(inv_r2, n) * inv_r;
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
            int n1, n2, n3;
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            Nabla_R[i] = Nabla_1_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe, cnt, coef);
        }
	}
    //Calculate the other elements
	fill_traceless_tensor(Nabla_R);
}


//Calculate the tensor Nabla 1/r for a give r(x,y,z) using recursive relation
void Nabla_r_traceless(double x, double y, double z, double * coef, double *Nabla_R)
{
    double r_2 = x*x + y*y + z*z;
    double inv_r2 = 1.0/r_2;
    double inv_r = sqrt(inv_r2);

    pow_x[0] = 1;
    pow_y[0] = 1;
    pow_z[0] = 1;
    pow_r2[0] = 1;

    for(int i=1; i<n_Max_rank+1; ++i){
        pow_x[i] = x*pow_x[i-1];
        pow_y[i] = y*pow_y[i-1];
        pow_z[i] = z*pow_z[i-1];
        pow_r2[i] = r_2*pow_r2[i-1];
    }

    //Calculate the first 2n+1 elements for each rank n
    int cnt = 0;
    int rank_dirct_calc = 2;
    if(n_Max_rank<rank_dirct_calc) rank_dirct_calc = n_Max_rank;
    for(int n = 0; n<rank_dirct_calc+1; ++n) {
        //Use the formula to calculate nabla_1_over_r
        double r_coe = order_minus_one[n] * pow(inv_r2, n) * inv_r;
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
            int n1, n2, n3;
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            Nabla_R[i] = Nabla_1_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe, cnt, coef);
        }
    }
    for(int n = rank_dirct_calc+1; n<n_Max_rank+1; ++n) {
        //Use the recursive relation to calculate nabla_1_over_r
        double inv_n = 1.0/n;
//        for (int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i) {
        for (int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i) {
            int n1, n2, n3;
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];
            double T1 = 0;  //First term in the recursive relation
            double T2 = 0;  //Second term in the recursive relation

            if (n1>0) {
                int idx = find_index[n1-1][n2][n3];
                T1 += x*Nabla_R[idx]*inv_Factorial[n1-1]*inv_Factorial[n2]*inv_Factorial[n3];
            }
            if (n2>0) {
                int idx = find_index[n1][n2-1][n3];
                T1 += y*Nabla_R[idx]*inv_Factorial[n1]*inv_Factorial[n2-1]*inv_Factorial[n3];;
            }
            if (n3>0) {
                int idx = find_index[n1][n2][n3-1];
                T1 += z*Nabla_R[idx]*inv_Factorial[n1]*inv_Factorial[n2]*inv_Factorial[n3-1];;
            }
            T1 *= (2*n-1);
            if (n1>1) {
                int idx = find_index[n1-2][n2][n3];
                T2 += Nabla_R[idx]*inv_Factorial[n1-2]*inv_Factorial[n2]*inv_Factorial[n3];;
            }
            if (n2>1) {
                int idx = find_index[n1][n2-2][n3];
                T2 += Nabla_R[idx]*inv_Factorial[n1]*inv_Factorial[n2-2]*inv_Factorial[n3];;
            }
            //n3<2 for independent elements of a traceless totally symmetric tensor
            T2 *= (n-1);
            //Recursive relation
            Nabla_R[i] = -1*(T1+T2)*inv_n*inv_r2*Factorial[n1]*Factorial[n2]*Factorial[n3];
        }
    }
    //Calculate the other elements
	fill_traceless_tensor(Nabla_R);
}


//Calculate the trace, used in function detracer.
double tracer(double * tensor, int n1, int n2, int n3, int m){
    double trace = 0;
    for(int k1=0; k1<m+1; ++k1){
        for(int k2=0; k2<m-k1+1; ++k2){
            int k3 = m-k1-k2;
            int idx = find_index[n1+2*k1][n2+2*k2][n3+2*k3];
            trace += tensor[idx]*Factorial[m]*inv_Factorial[k1]*inv_Factorial[k2]*inv_Factorial[k3];
        }
    }
    return trace;
}

//Detracer operator, convert a totally symmetric tensor into a traceless totally symmetric tensor, NOT used in the following code.
void detracer(double * symmetric_tensor, double * traceless_tensor){
    for(int n=0; n<n_Max_rank+1; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
            int n1, n2, n3;
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];
            traceless_tensor[i] = 0;
            for(int m1=0; m1<=(n1/2);++m1){
                for(int m2=0; m2<=(n2/2); ++m2){
                    for(int m3=0; m3<=(n3/2);++m3){
                        int m = m1+m2+m3;
                        int idx1, idx2, idx3;
                        idx1 = n1-2*m1;
                        idx2 = n2-2*m2;
                        idx3 = n3-2*m3;

                        traceless_tensor[i] += order_minus_one[m]*Factorial_odd[n-m]
                          *combination_HH[n1][m1]*combination_HH[n2][m2]*combination_HH[n3][m3]
                          *tracer(symmetric_tensor, idx1, idx2, idx3, m);

                    }
                }
            }
        }
    }
    fill_traceless_tensor(traceless_tensor);
}

//Detracer operator, convert a totally symmetric tensor into a traceless totally symmetric tensor, NOT used in the following code.
void detracer_direct(double * symmetric_tensor, double * traceless_tensor){
    for(int i=0; i<Number_of_total_element;++i){
        int n1, n2, n3, n;
		n1 = index_n1[i];
		n2 = index_n2[i];
		n3 = index_n3[i];
		n = n1 + n2 + n3;
		traceless_tensor[i] = 0;
		for(int m1=0; m1<=(n1/2);++m1){
            for(int m2=0; m2<=(n2/2); ++m2){
                for(int m3=0; m3<=(n3/2);++m3){
                    int m = m1+m2+m3;
                    int idx1, idx2, idx3;
                    idx1 = n1-2*m1;
                    idx2 = n2-2*m2;
                    idx3 = n3-2*m3;
                    traceless_tensor[i] += order_minus_one[m]*Factorial_odd[n-m]
                          *combination_HH[n1][m1]*combination_HH[n2][m2]*combination_HH[n3][m3]
                          *tracer(symmetric_tensor, idx1, idx2, idx3, m);

                }
            }
		}
    }
}

//contraction of two tensors of the same rank
double Contraction_equal_rank(double * Tensor1, double * Tensor2, int n){
    double cntr = 0;
    for(int i = n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
        cntr += combination_coef[i] * Tensor1[i] * Tensor2[i];
    }
    return cntr;
}

//Contraction of two tensors. The high rank one is traceless totally symmetric, the low rank one is totally symmetric, the result is traceless totally symmetric
void Contraction_traceless(double *High_rank_Tensor, double *Low_rank_Tensor, double *HL_rank_Tensor, int m, int n)
{
	// m: High rank; n: Low rank; k: Final rank, k=m-n
	int k = m - n;

	memset(HL_rank_Tensor, 0, Number_of_total_element*sizeof(double));

	//	started index and End index
	int Started_index_HL = n_Rank_Multipole_Start_Position[k];
	int End_index_HL = n_Rank_Multipole_Start_Position[k + 1];

	int Started_index_L = n_Rank_Multipole_Start_Position[n];
	int End_index_L = n_Rank_Multipole_Start_Position[n + 1];

	for (int j = Started_index_HL; j < Started_index_HL + 2 * k + 1; j++)
	{
	    int k1, k2, k3;
		k1 = index_n1[j];
		k2 = index_n2[j];
		k3 = index_n3[j];

		for (int i = Started_index_L; i < End_index_L; i++)
		{
		    int n1, n2, n3;
			n1 = index_n1[i];
			n2 = index_n2[i];
			n3 = index_n3[i];

            //find index by (m1,m2,m3)
			int index = find_index[k1+n1][k2+n2][k3+n3];

			HL_rank_Tensor[j] += combination_coef[i] * High_rank_Tensor[index] * Low_rank_Tensor[i];
		}
	}

	for (int j = Started_index_HL + 2 * k + 1; j < End_index_HL; j++)
	{
	    int k1, k2, k3;

		k1 = index_n1[j];
		k2 = index_n2[j];
		k3 = index_n3[j];

		int index_a = find_index[k1 + 2][k2][k3 - 2];
		int index_b = find_index[k1][k2 + 2][k3 - 2];

		HL_rank_Tensor[j] = -HL_rank_Tensor[index_a] - HL_rank_Tensor[index_b];
	}

}

//Take derivative before contraction of two tensors of the same rank
//void Contraction_dr(double * Tensor1, double * Tensor2, int n, double x, double y, double z, double &cx, double &cy, double &cz){
//    double tmp[3] = {0,0,0};
//    double coef[3] = {1/x, 1/y, 1/z};
//    int ni[3];
//    for(int i = n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
//        ni[0] = index_n1[i];
//        ni[1] = index_n2[i];
//        ni[2] = index_n3[i];
//
////        for(int id=0; id<3; ++id) {
////                if (ni[id]>0) tmp[id] += coef[id]*ni[id]*combination_coef[i] * Tensor1[i] * Tensor2[i];
////        }
//        double comb = combination_coef[i] * Tensor1[i] * Tensor2[i];
//        for(int id=0; id<3; ++id) tmp[id] += coef[id]*ni[id]*comb;
//    }
//
//    cx += tmp[0];
//    cy += tmp[1];
//    cz += tmp[2];
//}

void Contraction_dr(double * Tensor1, double * Tensor2, int n, double x, double y, double z, double &cx, double &cy, double &cz){
    double tmp[3] = {0,0,0};
    double coef[3] = {1/x, 1/y, 1/z};
    int ni[3];
    for(int i = 1; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
        ni[0] = index_n1[i];
        ni[1] = index_n2[i];
        ni[2] = index_n3[i];

//        for(int id=0; id<3; ++id) {
//                if (ni[id]>0) tmp[id] += coef[id]*ni[id]*combination_coef[i] * Tensor1[i] * Tensor2[i];
//        }
        double comb = combination_coef[i] * Tensor1[i] * Tensor2[i];
        for(int id=0; id<3; ++id) tmp[id] += ni[id]*coef[id]*comb;
    }

    cx += tmp[0];
    cy += tmp[1];
    cz += tmp[2];

    if (x==0) cx = 0;
    if (y==0) cy = 0;
    if (z==0) cz = 0;
}


void Contraction(double *High_rank_Tensor, double *Low_rank_Tensor, double *HL_rank_Tensor, int m, int n)
{
	// m: High rank; n: Low rank; k: Final rank, k=m-n
	int k = m - n;
	int m1, m2, m3, n1, n2, n3, k1, k2, k3;
	int i, j;

	memset(HL_rank_Tensor, 0, Number_of_total_element*sizeof(double));

	//	started index and End index
	int Started_index_HL = n_Rank_Multipole_Start_Position[k];
	int End_index_HL = n_Rank_Multipole_Start_Position[k + 1];

	int Started_index_L = n_Rank_Multipole_Start_Position[n];
	int End_index_L = n_Rank_Multipole_Start_Position[n + 1];

	for (j = Started_index_HL; j < End_index_HL; j++)
	{
		k1 = index_n1[j];
		k2 = index_n2[j];
		k3 = index_n3[j];

		for (i = Started_index_L; i < End_index_L; i++)
		{
			n1 = index_n1[i];
			n2 = index_n2[i];
			n3 = index_n3[i];

			m1 = k1 + n1;
			m2 = k2 + n2;
			m3 = k3 + n3;

			int index = find_index[m1][m2][m3];

			HL_rank_Tensor[j] += Factorial[n] / (Factorial[n1] * Factorial[n2] * Factorial[n3]) * High_rank_Tensor[index] * Low_rank_Tensor[i];
		}
	}
}
