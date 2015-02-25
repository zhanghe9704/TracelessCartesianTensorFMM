/**********************************
FMMKernel.cpp
Define the Kernel functions for the multiple level fast multipole algorithm using tensors

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/



#include "head.hpp"


//This function is used in the "Charge_to_Multipole" function.
void Charge_2_M_pre(double Charge, double Multipole_x, double Multipole_y, double Multipole_z, double Charge_x, double Charge_y, double Charge_z, double *Born_Multipole_pre)
{
	double x = Charge_x - Multipole_x;
	double y = Charge_y - Multipole_y;
	double z = Charge_z - Multipole_z;

    Symmetric_Tensor(x,y,z,Born_Multipole_pre);
}

void Charge_to_Multipole(Box & box, double *Charge, double *Charge_x, double *Charge_y, double *Charge_z, unsigned long int *ptclist, double *Born_Multipole){

	double * temp_M = scratch;
	memset(temp_M, 0, Number_of_total_element*sizeof(double));
	memset(Born_Multipole, 0, Number_of_total_element*sizeof(double));

    unsigned long int idx = box.first_ptcl;
    double Multipole_x = box.center[0];
    double Multipole_y = box.center[1];
    double Multipole_z = box.center[2];
	for (unsigned int i = 0; i < box.n_ptcl; ++i)
	{
		Charge_2_M_pre(Charge[idx], Multipole_x, Multipole_y, Multipole_z, Charge_x[idx], Charge_y[idx], Charge_z[idx], temp_M);
		for (int j = 0; j < Number_of_total_element; ++j)
		{
			Born_Multipole[j] += temp_M[j];
		}
		idx = ptclist[idx];
	}

	for(int n=0; n<n_Max_rank+1; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1];++i) Born_Multipole[i] *= inv_Fact_order_minus_one[n];
	}
}


void Multipole_to_Multipole(double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double * multipole_coef, double *Old_M, double *New_M)
{
	int cnt = 0;

	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n1 = index_n1[i];
		int n2 = index_n2[i];
		int n3 = index_n3[i];

		New_M[i] = 0;
		for (int m1 = 0; m1 <= n1; m1++)
		{
			for (int m2 = 0; m2 <= n2; m2++)
			{
				for (int m3 = 0; m3 <= n3; m3++)
				{
				    int coef = 1;
					if (new_x<old_x) coef *= (1 - (m1 % 2) * 2);
                    if (new_y<old_y) coef *= (1 - (m2 % 2) * 2);
                    if (new_z<old_z) coef *= (1 - (m3 % 2) * 2);

					New_M[i] += multipole_coef[cnt]*coef*Old_M[Find_index(n1 - m1, n2 - m2, n3 - m3)];
					++cnt;
				}
			}
		}
	}
}

void multipole_to_multipole_coef(double boxsize, double * multipole_coef){
    int cnt = 0;
    for (int i = 0; i < Number_of_total_element; i++)
	{
		int n1 = index_n1[i];
		int n2 = index_n2[i];
		int n3 = index_n3[i];

		int n = n1 + n2 + n3;

		for (int m1 = 0; m1 <= n1; m1++)
		{
			for (int m2 = 0; m2 <= n2; m2++)
			{
				for (int m3 = 0; m3 <= n3; m3++)
				{
					int m = m1 + m2 + m3;
					multipole_coef[cnt] = Factorial[n - m] * inv_Factorial [n] * combination[n1][m1] * combination[n2][m2] * combination[n3][m3] * pow(0.5*boxsize, m);
					++cnt;
				}
			}
		}
	}
}

void update_multipole_to_multipole_coef(double * multipole_coef){
    int cnt = 0;
    for (int i = 0; i < Number_of_total_element; i++){
		int n1 = index_n1[i];
		int n2 = index_n2[i];
		int n3 = index_n3[i];

		for (int m1 = 0; m1 <= n1; m1++)
		{
			for (int m2 = 0; m2 <= n2; m2++)
			{
				for (int m3 = 0; m3 <= n3; m3++)
				{
					int m = m1 + m2 + m3;
					multipole_coef[cnt] *= pow(2,m);
					++cnt;
				}
			}
		}
	}
};

void Calc_Nabla_R(double boxsize, double * Nabla_R){
    int shift;
    for(int i3=2; i3<4; ++i3){
        for(int i2=0; i2<i3+1; ++i2){
            for(int i1=0; i1<i2+1; ++i1){
                shift = nabla_idx[i1][i2][i3-2];
                Nabla_r_traceless(i1*boxsize, i2*boxsize, i3*boxsize, Nabla_1_element_r_coef, &Nabla_R[shift*Number_of_total_element]);
            }
        }
    }
}

void update_Nabla_R(double * Nabla_R){
    for(int n=0; n<n_Max_rank+1; ++n){
        double k = pow(2,n+1);
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
                for(int j=0; j<16; ++j)   Nabla_R[i+j*Number_of_total_element] *= k;
        }
    }
}

void Multipole_to_Local(double *Multipole_for_trans, double Multi_x, double Multi_y, double Multi_z, double Local_x, double Local_y, double Local_z, double * Saved_Nabla_R, double boxsize, double *M2L_translation)
{
	int idx[3], seq[3];
	int sign[3] = {1,1,1};

	if (Local_x<Multi_x) sign[0] = -1;
	if (Local_y<Multi_y) sign[1] = -1;
	if (Local_z<Multi_z) sign[2] = -1;

	idx[0] = round((Local_x-Multi_x)*sign[0]/boxsize);
	idx[1] = round((Local_y-Multi_y)*sign[1]/boxsize);
	idx[2] = round((Local_z-Multi_z)*sign[2]/boxsize);

	int change = sequence3(idx,seq);
	int shift = nabla_idx[idx[seq[0]]][idx[seq[1]]][idx[seq[2]]-2];

    double *Nabla_R = scratch2;

    for(int n=0; n<=n_Max_rank;++n){
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
            int coef, index, ni[3];
            ni[0] = index_n1[i];
            ni[1] = index_n2[i];
            ni[2] = index_n3[i];
            if (change==0) index = i;
            else index = Find_index(ni[seq[0]],ni[seq[1]], ni[seq[2]]);

            coef = 1;
            if (ni[0]%2==1)    coef *= sign[0];
            if (ni[1]%2==1)    coef *= sign[1];
            if (ni[2]%2==1)    coef *= sign[2];
            Nabla_R[i] = coef*Saved_Nabla_R[shift*Number_of_total_element+index];
        }
    }

    fill_traceless_tensor(Nabla_R);

	memset(M2L_translation, 0, Number_of_total_element*sizeof(double));

	double *M2L_temp = scratch;
	memset(M2L_temp, 0, Number_of_total_element*sizeof(double));
	for (int k = 0; k <= n_Max_rank; k++)
	{
		for (int l = 0; l <= n_Max_rank - k; l++)
		{
			Contraction_traceless(Nabla_R, Multipole_for_trans, M2L_temp, l + k, l);
            for (int i = n_Rank_Multipole_Start_Position[k]; i < n_Rank_Multipole_Start_Position[k]+2*k+1; ++i)
			{
				M2L_translation[i] += M2L_temp[i];
			}
		}
	}

	for(int n=0; n<=n_Max_rank; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
            M2L_translation[i] = M2L_translation[i]*inv_Factorial[n];
        }
	}

	fill_traceless_tensor(M2L_translation);
}

void Calc_Rho_Tensor(double boxsize, double * Rho_Tensor){
    Symmetric_Tensor(0.25*boxsize, 0.25*boxsize, 0.25*boxsize, Rho_Tensor);
    for(int ix=0; ix<2; ++ix){
        for (int iy = 0; iy<2; ++iy){
            for(int iz=0; iz<2; ++iz){
                int shift = 0;
                if(ix>0) shift += 4;
                if(iy>0) shift += 2;
                if(iz>0) shift += 1;

                if(shift>0){
                    for(int n=0; n<n_Max_rank+1; ++n){
                        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
                            int coef, n1, n2, n3;
                            n1 = index_n1[i];
                            n2 = index_n2[i];
                            n3 = index_n3[i];
                            coef = 1;
                            if (ix>0)    coef *= order_minus_one[n1];
                            if (iy>0)    coef *= order_minus_one[n2];
                            if (iz>0)    coef *= order_minus_one[n3];

                            Rho_Tensor[i+shift*Number_of_total_element] = coef*Rho_Tensor[i];
                        }
                    }
                }
            }
        }
    }
}

void update_Rho_Tensor(double * Rho_Tensor){
    for(int n=0; n<n_Max_rank+1; ++n){
        double k = pow(0.5,n);
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
                for(int j=0; j<8; ++j)   Rho_Tensor[i+j*Number_of_total_element] *= k;
        }
    }
}

void Local_to_Local(double *Old_Local, double Old_x, double Old_y, double Old_z, double New_x, double New_y, double New_z, double * Origin_Rho_Tensor, double *New_Local)
{
	memset(New_Local, 0, Number_of_total_element * sizeof(double));

	double *Rho_Tensor = scratch;

    int shift = 0;
    if (New_x<Old_x)    shift += 4;
    if (New_y<Old_y)    shift += 2;
    if (New_z<Old_z)    shift += 1;
    memcpy ( Rho_Tensor, &Origin_Rho_Tensor[shift*Number_of_total_element], Number_of_total_element*sizeof(double) );

	double *L2L_temp = scratch2;
	for (int n = 0; n <= n_Max_rank; n++)
	{
		for (int m = n; m <= n_Max_rank; m++)
		{
			memset(L2L_temp, 0, Number_of_total_element * sizeof(double));
			Contraction_traceless(Old_Local, Rho_Tensor, L2L_temp, m, m - n);
			for (int i = n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
                New_Local[i] += L2L_temp[i] * combination[m][m - n];
			}
		}
	}

	fill_traceless_tensor(New_Local);
}

double MultipolePotential(double *Multipole, double Multi_x, double Multi_y, double Multi_z, double Poten_x, double Poten_y, double Poten_z)
{
	double multipole_potential = 0;

	double *Nabla_R = scratch;
	memset(Nabla_R, 0, Number_of_total_element * sizeof(double));
	Nabla_r_traceless(Poten_x-Multi_x, Poten_y-Multi_y, Poten_z-Multi_z, Nabla_1_element_r_coef, Nabla_R);

	for (int n_rank = 0; n_rank <= n_Max_rank; n_rank++) multipole_potential += Contraction_equal_rank(Multipole, Nabla_R, n_rank);

	return multipole_potential;
}


//This function is used in the "Charge_to_Local" function.
void Charge_to_Local_pre(double q, double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double *L_expansion)
{
	double x = old_x - new_x;
	double y = old_y - new_y;
	double z = old_z - new_z;
    double r_2 = x * x + y * y + z * z;
    double inv_r2 = 1.0/r_2;
    double inv_r = sqrt(inv_r2);
	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n1 = index_n1[i];
		int n2 = index_n2[i];
		int n3 = index_n3[i];

		int n = n1 + n2 + n3;

		double r_coe = order_minus_one[n]* pow(inv_r2, n)*inv_r;

		L_expansion[i] = q * Nabla_1_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe);
	}
}

void Charge_to_Local_pre_traceless(double q, double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double *L_expansion)
{
	double x = old_x - new_x;
	double y = old_y - new_y;
	double z = old_z - new_z;
	double r_2 = x * x + y * y + z * z;
    double inv_r2 = 1.0/r_2;
    double inv_r = sqrt(inv_r2);

    int cnt=0;
	for(int n=0; n<n_Max_rank+1; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
            int n1, n2, n3;
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            double r_coe = order_minus_one[n] * pow(inv_r2, n) * inv_r;

            L_expansion[i] = q * Nabla_1_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe, cnt, Nabla_1_element_r_coef);
        }
	}
	fill_traceless_tensor(L_expansion);
}

void Charge_to_Local_traceless(Box &box, unsigned long int * ptclist, double *q, double *old_x, double *old_y, double *old_z, double new_x, double new_y, double new_z, double *L_expansion)
{
	double *temp_L = scratch;
	memset(temp_L, 0, Number_of_total_element*sizeof(double));
	memset(L_expansion, 0, Number_of_total_element*sizeof(double));

    unsigned long int idx = box.first_ptcl;

	for (unsigned int i = 0; i < box.n_ptcl; ++i)
	{
		Charge_to_Local_pre_traceless(q[idx], old_x[idx], old_y[idx], old_z[idx], new_x, new_y, new_z, temp_L);

        for(int n=0; n<n_Max_rank+1; ++n){
            for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
                L_expansion[i] += temp_L[i];
            }
        }
        idx = ptclist[idx];
	}

    for(int n=0; n<n_Max_rank+1; ++n){
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i) L_expansion[i] *= inv_Fact_order_minus_one[n];
    }

    fill_traceless_tensor(L_expansion);

}


double LocalPotential(double *Local_expan, double Local_x, double Local_y, double Local_z, double observer_x, double observer_y, double observer_z)
{
	double x = observer_x - Local_x;
	double y = observer_y - Local_y;
	double z = observer_z - Local_z;

	double L_Potential = 0;

	double *SymmeticTensor = scratch;
	memset(SymmeticTensor, 0, Number_of_total_element*sizeof(double));
	Symmetric_Tensor(x, y, z, SymmeticTensor);

	for (int i = 0; i <= n_Max_rank; i++) L_Potential += Contraction_equal_rank(Local_expan, SymmeticTensor,i);

	return L_Potential;
}
