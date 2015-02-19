/**********************************
Function.cpp
Define the math functions used in the kernel functions for the multiple level fast multipole algorithm using tensors

version 1.0
By He Huang & He Zhang, 12/29/2014

***********************************/


#include "head.hpp"

//int Find_index(int n1, int n2, int n3)
//{
//	int n = n1 + n2 + n3;
//	int index = n2 + (2 * n + 3 - n3)* n3 / 2;  // c++ index starting from 0
//
//	index = index + n_Rank_Multipole_Start_Position[n];
//	return index;
//}

int Find_index(int n1, int n2, int n3){return find_index[n1][n2][n3]+n_Rank_Multipole_Start_Position[n1 + n2 + n3];}

double combination(int n, int m)
{
	return Factorial[n] / (Factorial[m] * Factorial[n - m]);
}

double combination_HH(int n, int m)
{
	return Factorial[n] / (pow(2.0, m)* Factorial[m] * Factorial[n - 2 * m]);
}

void Symmetric_Tensor(double x, double y, double z, double *SymmeticTensor)
{
	for (int i = 0; i < Number_of_total_element; i++)
	{
		SymmeticTensor[i] =
			pow(x, index_n1[i]) * pow(y, index_n2[i]) * pow(z, index_n3[i]);
	}
}

double Nabla_1_element_r(int n1, int n2, int n3, int n, double x, double y, double z, double r_2, double r_coe)
{
	// reason of inputting "n, r^2, r_coe" is to accelerate calculation for the other fuctions
	// r_2 = x^2 + y^2 + z^2
	// r_coe = (-1)^n * r^(-2n)
	// n = n1+n2+n3

	double one_element = 0;
	int m1, m2, m3, m;

	for (m1 = 0; m1 <= (n1 / 2); m1++)
	{
		for (m2 = 0; m2 <= (n2 / 2); m2++)
		{
			for (m3 = 0; m3 <= (n3 / 2); m3++)
			{
				m = m1 + m2 + m3;

				one_element = one_element +
					(1 - (m % 2) * 2) * combination_HH(n1, m1) * combination_HH(n2, m2) *combination_HH(n3, m3) *
					pow(r_2, m) * Factorial_odd[n - m] *
					pow(x, n1 - 2 * m1) * pow(y, n2 - 2 * m2) * pow(z, n3 - 2 * m3);
			}
		}
	}

	one_element = one_element * r_coe;

	return one_element;
}

void Nabla_r_direct(double end_x, double end_y, double end_z, double begin_x, double begin_y, double begin_z, double *Nabla_R)
{
	double x = end_x - begin_x;
	double y = end_y - begin_y;
	double z = end_z - begin_z;

	for (int i = 0; i < Number_of_total_element; i++)
	{
		int n1, n2, n3;
		n1 = index_n1[i];
		n2 = index_n2[i];
		n3 = index_n3[i];

		int n = n1 + n2 + n3;
		double r_2 = x*x + y*y + z*z;
		double r_coe = (1 - (n % 2) * 2) / pow(r_2, n) / sqrt(r_2);

		double temp = Nabla_1_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe);

		Nabla_R[i] = temp;
	}
}

void Nabla_r_traceless(double end_x, double end_y, double end_z, double begin_x, double begin_y, double begin_z, double *Nabla_R)
{
	double x = end_x - begin_x;
	double y = end_y - begin_y;
	double z = end_z - begin_z;

	for (int rank_n = 0; rank_n <= n_Max_rank; rank_n++)
	{
		for (int j = 0; j < (2 * rank_n + 1); j++)
		{
			int n1, n2, n3;
			n1 = index_n1[n_Rank_Multipole_Start_Position[rank_n] + j];
			n2 = index_n2[n_Rank_Multipole_Start_Position[rank_n] + j];
			n3 = index_n3[n_Rank_Multipole_Start_Position[rank_n] + j];

			double r_2 = x*x + y*y + z*z;
			double r_coe = (1 - (rank_n % 2) * 2) / pow(r_2, rank_n) / sqrt(r_2);

			double temp = Nabla_1_element_r(n1, n2, n3, rank_n, x, y, z, r_2, r_coe);

			Nabla_R[n_Rank_Multipole_Start_Position[rank_n] + j] = temp;
		}

		for (int j = 2 * rank_n + 1; j < (rank_n + 1)*(rank_n + 2) / 2; j++)
		{
			int n1, n2, n3;
			n1 = index_n1[j + n_Rank_Multipole_Start_Position[rank_n]];
			n2 = index_n2[j + n_Rank_Multipole_Start_Position[rank_n]];
			n3 = index_n3[j + n_Rank_Multipole_Start_Position[rank_n]];

			int index_a = Find_index(n1 + 2, n2, n3 - 2);
			int index_b = Find_index(n1, n2 + 2, n3 - 2);

			Nabla_R[j + n_Rank_Multipole_Start_Position[rank_n]] = -Nabla_R[index_a] - Nabla_R[index_b];
		}
	}
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

			int index = Find_index(m1, m2, m3);

			HL_rank_Tensor[j] += Factorial[n] / (Factorial[n1] * Factorial[n2] * Factorial[n3]) * High_rank_Tensor[index] * Low_rank_Tensor[i];
		}
	}
}

void Contraction_traceless(double *High_rank_Tensor, double *Low_rank_Tensor, double *HL_rank_Tensor, int m, int n)
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

	for (j = Started_index_HL; j < Started_index_HL + 2 * k + 1; j++)
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

			int index = Find_index(m1, m2, m3);

			HL_rank_Tensor[j] += Factorial[n] / (Factorial[n1] * Factorial[n2] * Factorial[n3]) * High_rank_Tensor[index] * Low_rank_Tensor[i];
		}
	}

	for (int j = Started_index_HL + 2 * k + 1; j < End_index_HL; j++)
	{
		k1 = index_n1[j];
		k2 = index_n2[j];
		k3 = index_n3[j];

		int index_a = Find_index(k1 + 2, k2, k3 - 2);
		int index_b = Find_index(k1, k2 + 2, k3 - 2);

		HL_rank_Tensor[j] = -HL_rank_Tensor[index_a] - HL_rank_Tensor[index_b];
	}

}


