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

	for (int i = 0; i < Number_of_total_element; i++)
	{
		Born_Multipole_pre[i] = Charge * pow(x, index_n1[i]) * pow(y, index_n2[i]) * pow(z, index_n3[i]);
	}
}



//Charge to multipole
void Charge_to_Multipole(int Charge_Number, double *Charge, double Multipole_x, double Multipole_y, double Multipole_z, double *Charge_x, double *Charge_y, double *Charge_z, double *Born_Multipole){

	double * temp_M = scratch;
	memset(temp_M, 0, Number_of_total_element*sizeof(double));
	memset(Born_Multipole, 0, Number_of_total_element*sizeof(double));

	for (int i = 0; i < Charge_Number; i++)
	{
		Charge_2_M_pre(Charge[i], Multipole_x, Multipole_y, Multipole_z, Charge_x[i], Charge_y[i], Charge_z[i], temp_M);
		for (int j = 0; j < Number_of_total_element; j++)
		{
			Born_Multipole[j] += temp_M[j];
		}

	}

	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n = index_n1[i] + index_n2[i] + index_n3[i];

		Born_Multipole[i] *= (1 - (n % 2) * 2) / Factorial[n];
	}

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

	for (int i = 0; i < Number_of_total_element; ++i)
	{
		int n = index_n1[i] + index_n2[i] + index_n3[i];

		Born_Multipole[i] *= (1 - (n % 2) * 2) / Factorial[n];
	}

}




void Multipole_to_Multipole(double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double *Old_M, double *New_M)
{
	double x = new_x - old_x;
	double y = new_y - old_y;
	double z = new_z - old_z;

	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n1 = index_n1[i];
		int n2 = index_n2[i];
		int n3 = index_n3[i];

		int n = n1 + n2 + n3;

		New_M[i] = 0;
		for (int m1 = 0; m1 <= n1; m1++)
		{
			for (int m2 = 0; m2 <= n2; m2++)
			{
				for (int m3 = 0; m3 <= n3; m3++)
				{
					int m = m1 + m2 + m3;

					double Power_of_xyz = pow(x, m1) * pow(y, m2) * pow(z, m3);

					int index_of_old_element = Find_index(n1 - m1, n2 - m2, n3 - m3);
					double Old_element = Old_M[index_of_old_element];

					New_M[i] += Factorial[n - m] / Factorial[n] * combination(n1, m1) * combination(n2, m2) * combination(n3, m3) * Power_of_xyz * Old_element;
				}
			}
		}
	}
}


void Multipole_to_Multipole(double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double * multipole_coef, double *Old_M, double *New_M)
{
	double coef;

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
				    coef = 1;
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
					multipole_coef[cnt] = Factorial[n - m] / Factorial[n] * combination(n1, m1) * combination(n2, m2) * combination(n3, m3) * pow(0.5*boxsize, m);
					++cnt;
				}
			}
		}
	}
}

void update_multipole_to_multipole_coef(double * multipole_coef){
    int cnt = 0;
    for (int i = 0; i < Number_of_total_element; i++)
	{
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



void Multipole_to_Local(double *Multipole_for_trans, double Multi_x, double Multi_y, double Multi_z, double Local_x, double Local_y, double Local_z, double *M2L_translation)
{
	//double x = Multi_x - Local_x;
	//double y = Multi_y - Local_y;
	//double z = Multi_z - Local_z;

	memset(M2L_translation, 0, Number_of_total_element*sizeof(double));

	double *Nabla_R = scratch2;
	memset(Nabla_R, 0, Number_of_total_element*sizeof(double));
	Nabla_r_traceless(Local_x, Local_y, Local_z, Multi_x, Multi_y, Multi_z, Nabla_R);
	//	Nabla.Nabla_r_traceless(Multi_x, Multi_y, Multi_z, Local_x, Local_y, Local_z, Nabla_R);

	double *M2L_temp = scratch;
	memset(M2L_temp, 0, Number_of_total_element*sizeof(double));
	for (int k = 0; k <= n_Max_rank; k++)
	{
		for (int l = 0; l <= n_Max_rank - k; l++)
		{
			Contraction_traceless(Nabla_R, Multipole_for_trans, M2L_temp, l + k, l);
			for (int i = n_Rank_Multipole_Start_Position[k]; i < n_Rank_Multipole_Start_Position[k + 1]; i++)
			{
				M2L_translation[i] += M2L_temp[i];
				//				M2L_translation[i] += M2L_temp[i] * (1 - ((k + l) % 2) * 2);
			}
		}
	}
	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n = index_n1[i] + index_n2[i] + index_n3[i];
		M2L_translation[i] = M2L_translation[i] / Factorial[n];
	}

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
                        int coef, n1, n2, n3;
                        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
                            n1 = index_n1[i];
                            n2 = index_n2[i];
                            n3 = index_n3[i];
                            coef = 1;
                            if (ix>0)    coef *= (1 - (n1 % 2) * 2);
                            if (iy>0)    coef *= (1 - (n2 % 2) * 2);
                            if (iz>0)    coef *= (1 - (n3 % 2) * 2);

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
        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
                for(int j=0; j<8; ++j)   Rho_Tensor[i+j*Number_of_total_element] *= pow(0.5,n);
        }
    }
}

void Local_to_Local(double *Old_Local, double Old_x, double Old_y, double Old_z, double New_x, double New_y, double New_z, double * Origin_Rho_Tensor, double *New_Local)
{
	memset(New_Local, 0, Number_of_total_element * sizeof(double));

	double *Rho_Tensor = scratch;
//	memset(Rho_Tensor, 0, Number_of_total_element * sizeof(double));
//	for(int n=0; n<n_Max_rank+1; ++n){
//        int coef, n1, n2, n3;
//        for(int i=n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n+1]; ++i){
//            n1 = index_n1[i];
//            n2 = index_n2[i];
//            n3 = index_n3[i];
//            coef = 1;
//            if (New_x<Old_x)    coef *= (1 - (n1 % 2) * 2);
//            if (New_y<Old_y)    coef *= (1 - (n2 % 2) * 2);
//            if (New_z<Old_z)    coef *= (1 - (n3 % 2) * 2);
//
//            Rho_Tensor[i] = coef*Origin_Rho_Tensor[i];
//        }
//    }

    int shift = 0;
    if (New_x<Old_x)    shift += 4;
    if (New_y<Old_y)    shift += 2;
    if (New_z<Old_z)    shift += 1;
    memcpy ( Rho_Tensor, &Origin_Rho_Tensor[shift*Number_of_total_element], Number_of_total_element*sizeof(double) );

//	double x = New_x - Old_x;
//	double y = New_y - Old_y;
//	double z = New_z - Old_z;
//
//	Symmetric_Tensor(x, y, z, Rho_Tensor);

	double *L2L_temp = scratch2;
	for (int n = 0; n <= n_Max_rank; n++)
	{
		for (int m = n; m <= n_Max_rank; m++)
		{
			double Factor = combination(m, m - n);
//			double *L2L_temp = scratch2;
			memset(L2L_temp, 0, Number_of_total_element * sizeof(double));
			Contraction_traceless(Old_Local, Rho_Tensor, L2L_temp, m, m - n);

//			for (int i = n_Rank_Multipole_Start_Position[n]; i < n_Rank_Multipole_Start_Position[n + 1]; i++)
//			{
//				New_Local[i] += L2L_temp[i] * Factor;
//			}

			for (int i = n_Rank_Multipole_Start_Position[n]; i<n_Rank_Multipole_Start_Position[n]+2*n+1; ++i){
                New_Local[i] += L2L_temp[i] * Factor;
			}
		}
	}

	int n1, n2, n3;
    for(int n=0; n<n_Max_rank+1; ++n){
        int Start_index = n_Rank_Multipole_Start_Position[n];
        int End_index = n_Rank_Multipole_Start_Position[n + 1];
        for(int i=Start_index+2*n+1; i<End_index; ++i){
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            int index_a = Find_index(n1 + 2, n2, n3 - 2);
            int index_b = Find_index(n1, n2 + 2, n3 - 2);

            New_Local[i] = -New_Local[index_a] - New_Local[index_b];
        }
    }


}

void Local_to_Local(double *Old_Local, double Old_x, double Old_y, double Old_z, double New_x, double New_y, double New_z, double *New_Local)
{
	memset(New_Local, 0, Number_of_total_element * sizeof(double));

	double *Rho_Tensor = scratch;
	memset(Rho_Tensor, 0, Number_of_total_element * sizeof(double));

	double x = New_x - Old_x;
	double y = New_y - Old_y;
	double z = New_z - Old_z;

	Symmetric_Tensor(x, y, z, Rho_Tensor);

	double *L2L_temp = scratch2;
	for (int n = 0; n <= n_Max_rank; n++)
	{
		for (int m = n; m <= n_Max_rank; m++)
		{
			double Factor = combination(m, m - n);
//			double *L2L_temp = scratch2;
			memset(L2L_temp, 0, Number_of_total_element * sizeof(double));
			Contraction_traceless(Old_Local, Rho_Tensor, L2L_temp, m, m - n);

			for (int i = n_Rank_Multipole_Start_Position[n]; i < n_Rank_Multipole_Start_Position[n + 1]; i++)
			{
				New_Local[i] += L2L_temp[i] * Factor;
			}
		}
	}

}

double MultipolePotential(double *Multipole, double Multi_x, double Multi_y, double Multi_z, double Poten_x, double Poten_y, double Poten_z)
{
	double multipole_potential = 0;

	double *Nabla_R = scratch;
	memset(Nabla_R, 0, Number_of_total_element * sizeof(double));
	double *CT_temp = scratch2;
	memset(CT_temp, 0, Number_of_total_element * sizeof(double));
	Nabla_r_traceless(Poten_x, Poten_y, Poten_z, Multi_x, Multi_y, Multi_z, Nabla_R);

	for (int n_rank = 0; n_rank <= n_Max_rank; n_rank++)
	{
		Contraction(Multipole, Nabla_R, CT_temp, n_rank, n_rank);
		multipole_potential += CT_temp[0];
	}

	return multipole_potential;
}


//This function is used in the "Charge_to_Local" function.
void Charge_to_Local_pre(double q, double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double *L_expansion)
{
	double x = old_x - new_x;
	double y = old_y - new_y;
	double z = old_z - new_z;

	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n1 = index_n1[i];
		int n2 = index_n2[i];
		int n3 = index_n3[i];

		int n = n1 + n2 + n3;
		double r_2 = x * x + y * y + z * z;
		double r_coe = (1 - (n % 2) * 2) / (pow(r_2, n) * sqrt(r_2));

		L_expansion[i] = q * Nabla_1_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe);
	}
}

void Charge_to_Local_pre_traceless(double q, double old_x, double old_y, double old_z, double new_x, double new_y, double new_z, double *L_expansion)
{
	double x = old_x - new_x;
	double y = old_y - new_y;
	double z = old_z - new_z;

    int n1, n2, n3;
	for(int n=0; n<n_Max_rank+1; ++n){
        int Start_index = n_Rank_Multipole_Start_Position[n];
        int End_index = n_Rank_Multipole_Start_Position[n + 1];
        for(int i=Start_index; i<Start_index+2*n+1; ++i){
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            int n = n1 + n2 + n3;
            double r_2 = x * x + y * y + z * z;
            double r_coe = (1 - (n % 2) * 2) / (pow(r_2, n) * sqrt(r_2));

            L_expansion[i] = q * Nabla_1_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe);
        }
        for(int i=Start_index+2*n+1; i<End_index; ++i){
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            int index_a = Find_index(n1 + 2, n2, n3 - 2);
            int index_b = Find_index(n1, n2 + 2, n3 - 2);

            L_expansion[i] = -L_expansion[index_a] - L_expansion[index_b];
        }

	}
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
            int Start_index = n_Rank_Multipole_Start_Position[n];
            for(int i=Start_index; i<Start_index+2*n+1; ++i){
                L_expansion[i] += temp_L[i];
            }
        }
        idx = ptclist[idx];
	}

    int n1, n2, n3;
    for(int n=0; n<n_Max_rank+1; ++n){
        int Start_index = n_Rank_Multipole_Start_Position[n];
        int End_index = n_Rank_Multipole_Start_Position[n + 1];
        for(int i=Start_index; i<Start_index+2*n+1; ++i){
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            int n = n1 + n2 + n3;
            L_expansion[i] *= (1 - (n % 2) * 2) / Factorial[n];
        }
        for(int i=Start_index+2*n+1; i<End_index; ++i){
            n1 = index_n1[i];
            n2 = index_n2[i];
            n3 = index_n3[i];

            int index_a = Find_index(n1 + 2, n2, n3 - 2);
            int index_b = Find_index(n1, n2 + 2, n3 - 2);

            L_expansion[i] = -L_expansion[index_a] - L_expansion[index_b];
        }
    }

}


void Charge_to_Local(Box &box, unsigned long int * ptclist, double *q, double *old_x, double *old_y, double *old_z, double new_x, double new_y, double new_z, double *L_expansion)
{
	double *temp_L = scratch;
	memset(temp_L, 0, Number_of_total_element*sizeof(double));
	memset(L_expansion, 0, Number_of_total_element*sizeof(double));

    unsigned long int idx = box.first_ptcl;
	for (unsigned int i = 0; i < box.n_ptcl; ++i)
	{
		Charge_to_Local_pre(q[idx], old_x[idx], old_y[idx], old_z[idx], new_x, new_y, new_z, temp_L);
		for (int i = 0; i < Number_of_total_element; i++)
		{
			L_expansion[i] += temp_L[i];
		}
        idx = ptclist[idx];
	}


	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n1 = index_n1[i];
		int n2 = index_n2[i];
		int n3 = index_n3[i];

		int n = n1 + n2 + n3;
		L_expansion[i] *= (1 - (n % 2) * 2) / Factorial[n];
	}
}


void Charge_to_Local(double charge_number, double *q, double *old_x, double *old_y, double *old_z, double new_x, double new_y, double new_z, double *L_expansion)
{
	double *temp_L = scratch;
	memset(temp_L, 0, Number_of_total_element*sizeof(double));
	memset(L_expansion, 0, Number_of_total_element*sizeof(double));

	for (int i = 0; i < charge_number; i++)
	{
		Charge_to_Local_pre(q[i], old_x[i], old_y[i], old_z[i], new_x, new_y, new_z, temp_L);
		for (int i = 0; i < Number_of_total_element; i++)
		{
			L_expansion[i] += temp_L[i];
		}
	}
	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n1 = index_n1[i];
		int n2 = index_n2[i];
		int n3 = index_n3[i];

		int n = n1 + n2 + n3;
		L_expansion[i] *= (1 - (n % 2) * 2) / Factorial[n];
	}
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

	double *Local_temp = scratch2;
	memset(Local_temp, 0, Number_of_total_element*sizeof(double));

	for (int i = 0; i <= n_Max_rank; i++)
	{
		Contraction_traceless(Local_expan, SymmeticTensor, Local_temp, i, i);
		L_Potential += Local_temp[0];
	}

	return L_Potential;
}
