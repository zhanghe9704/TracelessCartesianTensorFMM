#include "fmm_kernel.h"
#include "global.h"

//This function is used in the "charge_to_multipole" function.
void charge_to_multipole_pre(double charge, double multipole_x, double multipole_y, double multipole_z,
                             double charge_x, double charge_y, double charge_z, double *multipole_pre) {
	double x = charge_x - multipole_x;
	double y = charge_y - multipole_y;
	double z = charge_z - multipole_z;
    symmetric_tensor_c2m(charge, x, y, z, multipole_pre);
}

void charge_to_multipole(Box &box, double *charge, double *charge_x, double *charge_y, double *charge_z,
                         double *multipole){
	double * temp_multipole = g_scratch;
	memset(temp_multipole, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));
	memset(multipole, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));

    unsigned long int idx = box.first_ptcl;
    double multipole_x = box.center[0];
    double multipole_y = box.center[1];
    double multipole_z = box.center[2];
	for (unsigned int i = 0; i < box.n_ptcl; ++i) {
		charge_to_multipole_pre(charge[idx], multipole_x, multipole_y, multipole_z, charge_x[idx], charge_y[idx],
                                charge_z[idx], temp_multipole);
		for (int j = 0; j < TOTAL_ELEMENT_NUMBER; ++j) {
			multipole[j] += temp_multipole[j];
		}
		idx = g_ptclist[idx];
	}
	for(int n=0; n<MAX_RANK+1; ++n){
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n+1];++i) multipole[i] *= kPowerNegOneInvFactN[n];
	}
}

void multipole_to_multipole(double old_x, double old_y, double old_z, double new_x, double new_y, double new_z,
                            double *multipole_coef, double *old_mulitpole, double *new_multipole) {
	int cnt = 0;
	memset(new_multipole, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));
	for (int i = 0; i < TOTAL_ELEMENT_NUMBER; i++) {
		int n1 = kIndexN1[i];
		int n2 = kIndexN2[i];
		int n3 = kIndexN3[i];
		for (int m1 = 0; m1 <= n1; m1++) {
			for (int m2 = 0; m2 <= n2; m2++) {
				for (int m3 = 0; m3 <= n3; m3++) {
				    int coef = 1;
					if (new_x<old_x) coef *= (1 - (m1 % 2) * 2);
                    if (new_y<old_y) coef *= (1 - (m2 % 2) * 2);
                    if (new_z<old_z) coef *= (1 - (m3 % 2) * 2);
					new_multipole[i] += multipole_coef[cnt]*coef*old_mulitpole[kIndex[n1 - m1][n2 - m2][n3 - m3]];
					++cnt;
				}
			}
		}
	}
}

void multipole_to_multipole_coef(double boxsize, double * multipole_coef) {
    int cnt = 0;
    for (int i = 0; i < TOTAL_ELEMENT_NUMBER; i++) {
		int n1 = kIndexN1[i];
		int n2 = kIndexN2[i];
		int n3 = kIndexN3[i];
		int n = n1 + n2 + n3;
		for (int m1 = 0; m1 <= n1; m1++) {
			for (int m2 = 0; m2 <= n2; m2++) {
				for (int m3 = 0; m3 <= n3; m3++) {
					int m = m1 + m2 + m3;
					multipole_coef[cnt] = kFactorial[n - m] * kInvFactorial [n] * kCombination[n1][m1] *
                                          kCombination[n2][m2] * kCombination[n3][m3] * pow(0.5*boxsize, m);
					++cnt;
				}
			}
		}
	}
}

void update_multipole_to_multipole_coef(double * multipole_coef) {
    int cnt = 0;
    for (int i = 0; i < TOTAL_ELEMENT_NUMBER; i++) {
		int n1 = kIndexN1[i];
		int n2 = kIndexN2[i];
		int n3 = kIndexN3[i];
		for (int m1 = 0; m1 <= n1; m1++) {
			for (int m2 = 0; m2 <= n2; m2++) {
				for (int m3 = 0; m3 <= n3; m3++) {
					int m = m1 + m2 + m3;
					multipole_coef[cnt] *= kPowerTwo[m];
					++cnt;
				}
			}
		}
	}
}

void calc_nabla_r(double boxsize, double *nabla_r){
    int shift;
    for(int i3=2; i3<4; ++i3) {
        for(int i2=0; i2<i3+1; ++i2) {
            for(int i1=0; i1<i2+1; ++i1) {
                shift = kNablaIndex[i1][i2][i3-2];
                nabla_r_traceless(i1*boxsize, i2*boxsize, i3*boxsize, g_nabla_coef, &nabla_r[shift*TOTAL_ELEMENT_NUMBER]);
            }
        }
    }
}

void update_nabla_r(double *nabla_r){
    for(int n=0; n<MAX_RANK+1; ++n){
        double k = kPowerTwo[n+1];
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n+1]; ++i){
                for(int j=0; j<16; ++j)   nabla_r[i+j*TOTAL_ELEMENT_NUMBER] *= k;
        }
    }
}

void multipole_to_local(double *multipole, double multi_x, double multi_y, double multi_z, double local_x,
                        double local_y, double local_z, double *saved_nabla_r, double boxsize, double *local_expn) {
	int idx[3], seq[3];
	int sign[3] = {1,1,1};
	if (local_x<multi_x) sign[0] = -1;
	if (local_y<multi_y) sign[1] = -1;
	if (local_z<multi_z) sign[2] = -1;
	idx[0] = round((local_x-multi_x)*sign[0]/boxsize);
	idx[1] = round((local_y-multi_y)*sign[1]/boxsize);
	idx[2] = round((local_z-multi_z)*sign[2]/boxsize);
	int change = sequence3(idx,seq);
	int shift = kNablaIndex[idx[seq[0]]][idx[seq[1]]][idx[seq[2]]-2];

    double *nabla_r = g_scratch2;
    for(int n=0; n<=MAX_RANK;++n){
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i){
            int coef, index, ni[3];
            ni[0] = kIndexN1[i];
            ni[1] = kIndexN2[i];
            ni[2] = kIndexN3[i];
            if (change==0) index = i;
            else index = kIndex[ni[seq[0]]][ni[seq[1]]][ni[seq[2]]];
            coef = 1;
            if (ni[0]%2==1)    coef *= sign[0];
            if (ni[1]%2==1)    coef *= sign[1];
            if (ni[2]%2==1)    coef *= sign[2];
            nabla_r[i] = coef*saved_nabla_r[shift*TOTAL_ELEMENT_NUMBER+index];
        }
    }
    fill_traceless_tensor(nabla_r);

	memset(local_expn, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));
	double *local_temp = g_scratch;
	memset(local_temp, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));
	for (int k = 0; k <= MAX_RANK; k++) {
		for (int l = 0; l <= MAX_RANK - k; l++) {
			contraction_traceless(nabla_r, multipole, local_temp, l + k, l);
            for (int i = kRankNTensorStart[k]; i < kRankNTensorStart[k]+2*k+1; ++i) {
				local_expn[i] += local_temp[i];
			}
		}
	}
	for(int n=0; n<=MAX_RANK; ++n) {
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i) {
            local_expn[i] = local_expn[i]*kInvFactorial[n];
        }
	}
	fill_traceless_tensor(local_expn);
}

void calc_rho_tensor(double boxsize, double *rho_tensor){
    symmetric_tensor(0.25*boxsize, 0.25*boxsize, 0.25*boxsize, rho_tensor);
    for(int ix=0; ix<2; ++ix) {
        for (int iy = 0; iy<2; ++iy) {
            for(int iz=0; iz<2; ++iz) {
                int shift = 0;
                if(ix>0) shift += 4;
                if(iy>0) shift += 2;
                if(iz>0) shift += 1;
                if(shift>0) {
                    for(int n=0; n<MAX_RANK+1; ++n) {
                        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n+1]; ++i) {
                            int coef, n1, n2, n3;
                            n1 = kIndexN1[i];
                            n2 = kIndexN2[i];
                            n3 = kIndexN3[i];
                            coef = 1;
                            if (ix>0)    coef *= kPowerNegOne[n1];
                            if (iy>0)    coef *= kPowerNegOne[n2];
                            if (iz>0)    coef *= kPowerNegOne[n3];
                            rho_tensor[i+shift*TOTAL_ELEMENT_NUMBER] = coef*rho_tensor[i];
                        }
                    }
                }
            }
        }
    }
}


void update_rho_tensor(double *rho_tensor) {
    for(int n=0; n<MAX_RANK+1; ++n) {
        double k = kPowerHalf[n];
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n+1]; ++i) {
                for(int j=0; j<8; ++j)   rho_tensor[i+j*TOTAL_ELEMENT_NUMBER] *= k;
        }
    }
}

void local_to_local(double *old_local, double old_x, double old_y, double old_z, double new_x, double new_y,
                    double new_z, double *origin_rho_tensor, double *new_local) {
	memset(new_local, 0, TOTAL_ELEMENT_NUMBER * sizeof(double));
	double *rho_tensor = g_scratch;
    int shift = 0;
    if (new_x<old_x)    shift += 4;
    if (new_y<old_y)    shift += 2;
    if (new_z<old_z)    shift += 1;
    memcpy (rho_tensor, &origin_rho_tensor[shift*TOTAL_ELEMENT_NUMBER], TOTAL_ELEMENT_NUMBER*sizeof(double));

	double *local_temp = g_scratch2;
	for (int n = 0; n <= MAX_RANK; n++) {
		for (int m = n; m <= MAX_RANK; m++) {
			contraction_traceless(old_local, rho_tensor, local_temp, m, m - n);
			for (int i = kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i) {
                new_local[i] += local_temp[i] * kCombination[m][m - n];
			}
		}
	}
	fill_traceless_tensor(new_local);
}

int multipole_to_charge(double *multipole, double multi_x, double multi_y, double multi_z, double obsv_x, double obsv_y,
                        double obsv_z, double &phi, double &ex, double &ey, double &ez) {
    double *nabla_r = g_scratch;
    memset(nabla_r, 0, TOTAL_ELEMENT_NUMBER * sizeof(double));
    nabla_r_traceless(obsv_x-multi_x, obsv_y-multi_y, obsv_z-multi_z, g_nabla_coef, nabla_r);

    if(g_flag!=Flag::FIELD) {
        double potential = 0;
        for (int n_rank = 0; n_rank <= MAX_RANK; n_rank++)
            potential += contraction_equal_rank(multipole, nabla_r, n_rank);
        phi += potential;
    }
    if(g_flag!=Flag::POTENTIAL) {
        double *field = g_scratch2;
        for (int n_rank = 0; n_rank < MAX_RANK; ++n_rank){
            contraction_traceless(nabla_r,multipole,field,n_rank+1,n_rank);
            ex -= field[1];
            ey -= field[2];
            ez -= field[3];
        }
    }
    return 0;
}

//This function is used in the "charge_to_local_traceless" function.
void charge_to_local_pre_traceless(double q, double old_x, double old_y, double old_z, double new_x, double new_y,
                                   double new_z, double *local_expn) {
	double x = old_x - new_x;
	double y = old_y - new_y;
	double z = old_z - new_z;
	double r_2 = x * x + y * y + z * z;
    double inv_r2 = 1.0/r_2;
    double inv_r = sqrt(inv_r2);
    int cnt=0;

    g_pow_x[0] = 1;
    g_pow_y[0] = 1;
    g_pow_z[0] = 1;
    g_pow_r2[0] = 1;
    for(int i=1; i<MAX_RANK+1; ++i) {
        g_pow_x[i] = x*g_pow_x[i-1];
        g_pow_y[i] = y*g_pow_y[i-1];
        g_pow_z[i] = z*g_pow_z[i-1];
        g_pow_r2[i] = r_2*g_pow_r2[i-1];
    }
	for(int n=0; n<MAX_RANK+1; ++n) {
        double r_coe = kPowerNegOne[n] * pow(inv_r2, n) * inv_r;
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i){
            int n1, n2, n3;
            n1 = kIndexN1[i];
            n2 = kIndexN2[i];
            n3 = kIndexN3[i];
            local_expn[i] = q * nabla_element_r(n1, n2, n3, n, x, y, z, r_2, r_coe, cnt, g_nabla_coef);
        }
	}
	fill_traceless_tensor(local_expn);
}

void charge_to_local_traceless(Box &box, double *q, double *old_x, double *old_y, double *old_z, double new_x,
                               double new_y, double new_z, double *local_expn) {
	double *local_temp = g_scratch;
	memset(local_temp, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));
	memset(local_expn, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));

    unsigned long int idx = box.first_ptcl;
	for (unsigned int i = 0; i < box.n_ptcl; ++i) {
		charge_to_local_pre_traceless(q[idx], old_x[idx], old_y[idx], old_z[idx], new_x, new_y, new_z, local_temp);
        for(int n=0; n<MAX_RANK+1; ++n) {
            for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i){
                local_expn[i] += local_temp[i];
            }
        }
        idx = g_ptclist[idx];
	}
    for(int n=0; n<MAX_RANK+1; ++n) {
        for(int i=kRankNTensorStart[n]; i<kRankNTensorStart[n]+2*n+1; ++i) local_expn[i] *= kPowerNegOneInvFactN[n];
    }
    fill_traceless_tensor(local_expn);
}

int local_to_charge(double *local_expn, double local_x, double local_y, double local_z, double obsv_x, double obsv_y,
                    double obsv_z, double &phi, double &ex, double &ey, double &ez) {
    double dx = obsv_x - local_x;
    double dy = obsv_y - local_y;
    double dz = obsv_z - local_z;
    double *tensor_rn = g_scratch;
    memset(tensor_rn, 0, TOTAL_ELEMENT_NUMBER*sizeof(double));
	symmetric_tensor(dx, dy, dz, tensor_rn);
	if(g_flag!=Flag::FIELD) {
        double potential = 0;
        for (int i = 0; i <= MAX_RANK; i++) potential += contraction_equal_rank(local_expn, tensor_rn,i);
        phi += potential;
	}
	if(g_flag!=Flag::POTENTIAL) {
        double field_x = 0;
        double field_y = 0;
        double field_z = 0;
        contraction_dr(local_expn, tensor_rn, MAX_RANK, dx, dy, dz, field_x, field_y, field_z);
        ex -= field_x;
        ey -= field_y;
        ez -= field_z;
	}
	return 0;
}

