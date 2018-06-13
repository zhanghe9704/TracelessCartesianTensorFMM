#include "global.h"

//output i3>=i2>=i1, if sequence changed return a number greater than 0;
int sequence3 (int x[3], int idx[3]) {
    int change = 0;
    for (int i=0; i<3; ++i) idx[i] = i;
    if (x[1]>x[2])  {idx[2]=1; idx[1]=2; ++change; }
    if (x[0]>x[idx[2]]) {idx[0]=idx[2]; idx[2]=0; ++change;}
    if (x[idx[0]]>x[idx[1]]) { int tmp = idx[0]; idx[0] = idx[1]; idx[1] = tmp; ++change;}
    return change;
}

//fill the totally symmetric tensor r^n
void fill_symmetric_tensor_r(double r2, double *symmetric_tensor) {
    for(int n=0; n<MAX_RANK+1; ++n){
        for(int i=kRankNTensorStart[n]+2*n+1; i<kRankNTensorStart[n+1]; ++i){
            int n1, n2, n3, index_a, index_b, index;
            n1 = kIndexN1[i];
            n2 = kIndexN2[i];
            n3 = kIndexN3[i];
            index_a = kIndex[n1 + 2][n2][n3 - 2];
            index_b = kIndex[n1][n2 + 2][n3 - 2];
            index = kIndex[n1][n2][n3-2];
            symmetric_tensor[i] = symmetric_tensor[index]*r2-symmetric_tensor[index_a]-symmetric_tensor[index_b];
        }
    }
}

//Calculate the totally symmetric tensor r^n
void symmetric_tensor(double x, double y, double z, double *tensor) {
    tensor[0] = 1;
    if (MAX_RANK>0) {
        tensor[1] = x;
        tensor[2] = y;
        tensor[3] = z;
    }
    double r2 = x*x+y*y+z*z;
    for(int n=2; n<MAX_RANK+1; ++n) {
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i) {
            int n1 = kIndexN1[i];
            int n2 = kIndexN2[i];
            int n3 = kIndexN3[i];
            if (n1>0) {
                int index = kIndex[n1-1][n2][n3];
                tensor[i] = tensor[index]*x;
            }
            else {
                int index = kIndex[n1][n2-1][n3];
                tensor[i] = tensor[index]*y;
            }
        }
    }
    fill_symmetric_tensor_r(r2, tensor);
}

//Calculate the totally symmetric tensor r^n times charge Q for charge to multipole calculation
void symmetric_tensor_c2m(double q, double x, double y, double z, double *symmetric_tensor){
    symmetric_tensor[0] = q;
    if (MAX_RANK>0){
        symmetric_tensor[1] = x*q;
        symmetric_tensor[2] = y*q;
        symmetric_tensor[3] = z*q;
    }
    double r2 = x*x+y*y+z*z;
    for(int n=2; n<MAX_RANK+1; ++n){
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i){
            int n1 = kIndexN1[i];
            int n2 = kIndexN2[i];
            int n3 = kIndexN3[i];
            if (n1>0){
                int index = kIndex[n1-1][n2][n3];
                symmetric_tensor[i] = symmetric_tensor[index]*x;
            }
            else{
                int index = kIndex[n1][n2-1][n3];
                symmetric_tensor[i] = symmetric_tensor[index]*y;
            }
        }
    }
    fill_symmetric_tensor_r(r2, symmetric_tensor);
}

//Using the value of the first 2*n+1 elements to calculate the value of the other elements for a traceless totally
//symmetric tensor
void fill_traceless_tensor(double *tensor){
    for(int n=0; n<=MAX_RANK; ++n){
        for(int i=kRankNTensorStart[n]+2*n+1; i<kRankNTensorStart[n+1]; ++i){
            int n1, n2, n3, index_a, index_b;
            n1 = kIndexN1[i];
            n2 = kIndexN2[i];
            n3 = kIndexN3[i];
            index_a = kIndex[n1 + 2][n2][n3 - 2];
            index_b = kIndex[n1][n2 + 2][n3 - 2];
            tensor[i] = -tensor[index_a] - tensor[index_b];
        }
	}
}

//Calculate the value of a specific element of the tensor/operator Nabla 1/r
double nabla_element_r(int n1, int n2, int n3, int n, double x, double y, double z, double r_2, double r_coe, int &cnt,
                       double *coef) {
	// reason of inputting "n, r^2, r_coe" is to accelerate calculation for the other functions
	// r_2 = x^2 + y^2 + z^2
	// r_coe = (-1)^n * r^(-2n)
	// n = n1+n2+n3

	double one_element = 0;
	for (int m1 = 0; m1 <= (n1 / 2); m1++){
		for (int m2 = 0; m2 <= (n2 / 2); m2++){
			for (int m3 = 0; m3 <= (n3 / 2); m3++){
				int m = m1 + m2 + m3;
				one_element += coef[cnt] * g_pow_r2[m]*g_pow_x[n1 - 2 * m1]*g_pow_y[n2 - 2 * m2]*g_pow_z[n3 - 2 * m3];
				++cnt;
			}
		}
	}
	one_element *= r_coe;
	return one_element;
}

//Calculate the tensor Nabla 1/r for a give r(x,y,z) using recursive relation
void nabla_r_traceless(double x, double y, double z, double *coef, double *nabla_r)
{
    double r_2 = x*x + y*y + z*z;
    double inv_r2 = 1.0/r_2;
    double inv_r = sqrt(inv_r2);

    g_pow_x[0] = 1;
    g_pow_y[0] = 1;
    g_pow_z[0] = 1;
    g_pow_r2[0] = 1;
    for(int i=1; i<MAX_RANK+1; ++i){
        g_pow_x[i] = x*g_pow_x[i-1];
        g_pow_y[i] = y*g_pow_y[i-1];
        g_pow_z[i] = z*g_pow_z[i-1];
        g_pow_r2[i] = r_2*g_pow_r2[i-1];
    }
    //Calculate the first 2n+1 elements for each rank n
    int cnt = 0;
    int rank_dirct_calc = MAX_RANK;
//    int rank_dirct_calc = 2;
    if(MAX_RANK<rank_dirct_calc) rank_dirct_calc = MAX_RANK;
    for(int n = 0; n<rank_dirct_calc+1; ++n) {
        //Use the formula to calculate nabla_1_over_r
        double r_coe = kPowerNegOne[n] * pow(inv_r2, n) * inv_r;
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i){
            int n1, n2, n3;
            n1 = kIndexN1[i];
            n2 = kIndexN2[i];
            n3 = kIndexN3[i];
            nabla_r[i] = nabla_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe, cnt, coef);
        }
    }
    for(int n = rank_dirct_calc+1; n<MAX_RANK+1; ++n) {
        //Use the recursive relation to calculate nabla_1_over_r
        double inv_n = 1.0/n;
//        for (int i=kRankNTensorStart[n]; i<kRankNTensorStart[n+1]; ++i) {
        for (int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i) {
            int n1, n2, n3;
            n1 = kIndexN1[i];
            n2 = kIndexN2[i];
            n3 = kIndexN3[i];
            double T1 = 0;  //First term in the recursive relation
            double T2 = 0;  //Second term in the recursive relation
            if (n1>0) {
                int idx = kIndex[n1-1][n2][n3];
                T1 += x*nabla_r[idx]*kInvFactorial[n1-1]*kInvFactorial[n2]*kInvFactorial[n3];
            }
            if (n2>0) {
                int idx = kIndex[n1][n2-1][n3];
                T1 += y*nabla_r[idx]*kInvFactorial[n1]*kInvFactorial[n2-1]*kInvFactorial[n3];
            }
            if (n3>0) {
                int idx = kIndex[n1][n2][n3-1];
                T1 += z*nabla_r[idx]*kInvFactorial[n1]*kInvFactorial[n2]*kInvFactorial[n3-1];
            }
            T1 *= (2*n-1);
            if (n1>1) {
                int idx = kIndex[n1-2][n2][n3];
                T2 += nabla_r[idx]*kInvFactorial[n1-2]*kInvFactorial[n2]*kInvFactorial[n3];
            }
            if (n2>1) {
                int idx = kIndex[n1][n2-2][n3];
                T2 += nabla_r[idx]*kInvFactorial[n1]*kInvFactorial[n2-2]*kInvFactorial[n3];
            }
            //n3<2 for independent elements of a traceless totally symmetric tensor
            T2 *= (n-1);
            //Recursive relation
            nabla_r[i] = -1*(T1+T2)*inv_n*inv_r2*kFactorial[n1]*kFactorial[n2]*kFactorial[n3];
        }
    }
    //Calculate the other elements
	fill_traceless_tensor(nabla_r);
}

//contraction of two tensors of the same rank
double contraction_equal_rank(double *Tensor1, double *Tensor2, int n){
    double cntr = 0;
    for(int i = kRankNTensorStart[n]; i<kRankNTensorStart[n+1]; ++i){
        cntr += g_combination_coef[i] * Tensor1[i] * Tensor2[i];
    }
    return cntr;
}

//Contraction of two tensors. The high rank one is traceless totally symmetric, the low rank one is totally symmetric,
//the result is traceless totally symmetric
void contraction_traceless(double *high_rank_tensor, double *low_rank_tensor, double *output_tensor, int m, int n) {
	// m: High rank; n: Low rank; k: Final rank, k=m-n
	int k = m - n;
	memset(output_tensor, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));

	//	Start index and end index for the output tensor
	int output_start = kRankNTensorStart[k];
	int output_end = kRankNTensorStart[k + 1];
    //	Start index and end index for the low rank tensor
	int low_start = kRankNTensorStart[n];
	int low_end = kRankNTensorStart[n + 1];

	for (int j = output_start; j < output_start + 2 * k + 1; j++) {
	    int k1, k2, k3;
		k1 = kIndexN1[j];
		k2 = kIndexN2[j];
		k3 = kIndexN3[j];
		for (int i = low_start; i < low_end; i++) {
		    int n1, n2, n3;
			n1 = kIndexN1[i];
			n2 = kIndexN2[i];
			n3 = kIndexN3[i];
			int index = kIndex[k1+n1][k2+n2][k3+n3];
			output_tensor[j] += g_combination_coef[i] * high_rank_tensor[index] * low_rank_tensor[i];
		}
	}
	for (int j = output_start + 2 * k + 1; j < output_end; j++) {
	    int k1, k2, k3;
		k1 = kIndexN1[j];
		k2 = kIndexN2[j];
		k3 = kIndexN3[j];
		int index_a = kIndex[k1 + 2][k2][k3 - 2];
		int index_b = kIndex[k1][k2 + 2][k3 - 2];
		output_tensor[j] = -output_tensor[index_a] - output_tensor[index_b];
	}
}

void contraction_dr(double *tensor1, double *tensor2, int n, double x, double y, double z, double &cx, double &cy,
                    double &cz) {
    double tmp[3] = {0,0,0};
    double coef[3] = {1/x, 1/y, 1/z};
    int ni[3];
    for(int i = 1; i<kRankNTensorStart[n+1]; ++i){
        ni[0] = kIndexN1[i];
        ni[1] = kIndexN2[i];
        ni[2] = kIndexN3[i];
        double comb = g_combination_coef[i] * tensor1[i] * tensor2[i];
        for(int id=0; id<3; ++id) tmp[id] += ni[id]*coef[id]*comb;
    }
    cx += tmp[0];
    cy += tmp[1];
    cz += tmp[2];
    if (x==0) cx = 0;
    if (y==0) cy = 0;
    if (z==0) cz = 0;
}

//******************************************
// The following functions are not used.
//******************************************

//Contraction of two tensors of rank m and n. m>=n. The output tensor has rank m-n.
void contraction(double *high_rank_tensor, double *low_rank_tensor, double *output_tensor, int m, int n) {
	// m: High rank; n: Low rank; k: Final rank, k=m-n
	int k = m - n;
	int m1, m2, m3, n1, n2, n3, k1, k2, k3;
	int i, j;
	memset(output_tensor, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));

	//	Start index and end index for the output tensor
	int output_start = kRankNTensorStart[k];
	int output_end = kRankNTensorStart[k + 1];
    //	Start index and end index for the low rank tensor
	int low_start = kRankNTensorStart[n];
	int low_end = kRankNTensorStart[n + 1];
	for (j = output_start; j < output_end; j++) {
		k1 = kIndexN1[j];
		k2 = kIndexN2[j];
		k3 = kIndexN3[j];
		for (i = low_start; i < low_end; i++) {
			n1 = kIndexN1[i];
			n2 = kIndexN2[i];
			n3 = kIndexN3[i];

			m1 = k1 + n1;
			m2 = k2 + n2;
			m3 = k3 + n3;
			int index = kIndex[m1][m2][m3];
			output_tensor[j] += kFactorial[n] / (kFactorial[n1] * kFactorial[n2] * kFactorial[n3]) *
                                high_rank_tensor[index] * low_rank_tensor[i];
		}
	}
}

//Calculate the tensor Nabla 1/r for a give r(x,y,z)
void nabla_r_traceless_formula(double x, double y, double z, double *coef, double *nabla_r)
{
    double r_2 = x*x + y*y + z*z;
    double inv_r2 = 1.0/r_2;
    double inv_r = sqrt(inv_r2);

    g_pow_x[0] = 1;
    g_pow_y[0] = 1;
    g_pow_z[0] = 1;
    g_pow_r2[0] = 1;
    for(int i=1; i<MAX_RANK+1; ++i){
        g_pow_x[i] = x*g_pow_x[i-1];
        g_pow_y[i] = y*g_pow_y[i-1];
        g_pow_z[i] = z*g_pow_z[i-1];
        g_pow_r2[i] = r_2*g_pow_r2[i-1];
    }

    //Calculate the first 2n+1 elements for each rank n
    int cnt = 0;
	for(int n=0; n<MAX_RANK+1; ++n){
        double r_coe = kPowerNegOne[n] * pow(inv_r2, n) * inv_r;
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i){
            int n1, n2, n3;
            n1 = kIndexN1[i];
            n2 = kIndexN2[i];
            n3 = kIndexN3[i];
            nabla_r[i] = nabla_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe, cnt, coef);
        }
	}
    //Calculate the other elements
	fill_traceless_tensor(nabla_r);
}

//Calculate the trace, used in function detracer.
double tracer(double *tensor, int n1, int n2, int n3, int m){
    double trace = 0;
    for(int k1=0; k1<m+1; ++k1){
        for(int k2=0; k2<m-k1+1; ++k2){
            int k3 = m-k1-k2;
            int idx = kIndex[n1+2*k1][n2+2*k2][n3+2*k3];
            trace += tensor[idx]*kFactorial[m]*kInvFactorial[k1]*kInvFactorial[k2]*kInvFactorial[k3];
        }
    }
    return trace;
}

//Detracer operator, convert a totally symmetric tensor into a traceless totally symmetric tensor.
void detracer(double *symmetric_tensor, double *traceless_tensor){
    for(int n=0; n<MAX_RANK+1; ++n){
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i){
            int n1, n2, n3;
            n1 = kIndexN1[i];
            n2 = kIndexN2[i];
            n3 = kIndexN3[i];
            traceless_tensor[i] = 0;
            for(int m1=0; m1<=(n1/2);++m1){
                for(int m2=0; m2<=(n2/2); ++m2){
                    for(int m3=0; m3<=(n3/2);++m3){
                        int m = m1+m2+m3;
                        int idx1, idx2, idx3;
                        idx1 = n1-2*m1;
                        idx2 = n2-2*m2;
                        idx3 = n3-2*m3;
                        traceless_tensor[i] += kPowerNegOne[m]*kFactorialOdd[n-m]
                                              *kCombinationHH[n1][m1]*kCombinationHH[n2][m2]*kCombinationHH[n3][m3]
                                              *tracer(symmetric_tensor, idx1, idx2, idx3, m);
                    }
                }
            }
        }
    }
    fill_traceless_tensor(traceless_tensor);
}

//Detracer operator, convert a totally symmetric tensor into a traceless totally symmetric tensor.
void detracer_direct(double *symmetric_tensor, double *traceless_tensor){
    for(int i=0; i<TOTAL_ELEMENT_NUMBER;++i){
        int n1, n2, n3, n;
		n1 = kIndexN1[i];
		n2 = kIndexN2[i];
		n3 = kIndexN3[i];
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
                    traceless_tensor[i] += kPowerNegOne[m]*kFactorialOdd[n-m]
                                          *kCombinationHH[n1][m1]*kCombinationHH[n2][m2]*kCombinationHH[n3][m3]
                                          *tracer(symmetric_tensor, idx1, idx2, idx3, m);
                }
            }
		}
    }
}


