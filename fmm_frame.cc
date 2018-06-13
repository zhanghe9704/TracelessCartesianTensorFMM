#include <cassert>
#include "global.h"

//Calculate the multipole expansions for all boxes except for the root box
int calc_multipole(vector<Box> &tree, double *q, double *x, double *y, double *z) {
    double * multipole_coef = new double[kRankNTensorLength[MAX_RANK]];
    double current_childbox_size = tree[tree.size()-1].box_size;
    multipole_to_multipole_coef(current_childbox_size, multipole_coef);

    for(unsigned long int itr=tree.size()-1; itr>0; --itr) {
		if (tree[itr].n_child>0) {  //parent box, take summation of the multipole expansion of the child boxes
                for (int i=0; i<tree[itr].n_child; ++i) {
                        unsigned long int child_ptr = tree[itr].child[i];
                        //Translate the mltp of a child box to the parent box, use the first TOTAL_ELEMENT_NUMBER number
                        //of multipoles as g_scratch variable
                        if (tree[child_ptr].box_size>1.5*current_childbox_size) {
                            update_multipole_to_multipole_coef(multipole_coef);
                            current_childbox_size = tree[child_ptr].box_size;
                        }
                        multipole_to_multipole(tree[child_ptr].center[0],tree[child_ptr].center[1],
                                               tree[child_ptr].center[2],tree[itr].center[0],tree[itr].center[1],
                                               tree[itr].center[2], multipole_coef,
                                               &g_multipole_expns[child_ptr*TOTAL_ELEMENT_NUMBER], g_multipole_expns );
                        for (int j=0; j<TOTAL_ELEMENT_NUMBER; ++j) {
                            g_multipole_expns[itr*TOTAL_ELEMENT_NUMBER+j] += g_multipole_expns[j];
                        }
                }
		}
		else{	//childless box, calculate the multipole expansion from the charges
               charge_to_multipole(tree[itr], q, x, y, z, &g_multipole_expns[itr*TOTAL_ELEMENT_NUMBER]);
		}
	}
    delete[] multipole_coef;
	return 0;
}

int coulomb_potential(Box &box, double *x, double *y, double *z, double *q, double *phi,
                      double *ex, double *ey, double *ez) {
    unsigned long int obj_idx, src_idx;
    obj_idx = box.first_ptcl;
    while (obj_idx<TOTAL_PARTICLE_NUMBER) {
        src_idx = g_ptclist[obj_idx];
        while (src_idx<TOTAL_PARTICLE_NUMBER) {
            double dx = x[obj_idx]-x[src_idx];
            double dy = y[obj_idx]-y[src_idx];
            double dz = z[obj_idx]-z[src_idx];
            double r = 1/sqrt(dx*dx+dy*dy+dz*dz);
            if (g_flag!=Flag::FIELD) {
                phi[obj_idx] += q[src_idx]*r;
                phi[src_idx] += q[obj_idx]*r;
            }
            if (g_flag!=Flag::POTENTIAL) {
                double r3 = r*r*r;
                ex[obj_idx] += q[src_idx]*dx*r3;
                ey[obj_idx] += q[src_idx]*dy*r3;
                ez[obj_idx] += q[src_idx]*dz*r3;
                ex[src_idx] -= q[obj_idx]*dx*r3;
                ey[src_idx] -= q[obj_idx]*dy*r3;
                ez[src_idx] -= q[obj_idx]*dz*r3;
            }
            src_idx = g_ptclist[src_idx];
        }
        obj_idx = g_ptclist[obj_idx];
    }
    return 0;
}

int coulomb_potential(Box &obj_box, Box &src_box, double *x, double *y, double *z, double *q,
                      double *phi, double *ex, double *ey, double *ez) {
    unsigned long int obj_idx, src_idx;
    obj_idx = obj_box.first_ptcl;
    while (obj_idx<TOTAL_PARTICLE_NUMBER) {
        src_idx = src_box.first_ptcl;
        while (src_idx<TOTAL_PARTICLE_NUMBER) {
            double dx = x[obj_idx]-x[src_idx];
            double dy = y[obj_idx]-y[src_idx];
            double dz = z[obj_idx]-z[src_idx];
            double r = 1/sqrt(dx*dx+dy*dy+dz*dz);
            if (g_flag!=Flag::FIELD) {
                phi[obj_idx] += q[src_idx]*r;
                phi[src_idx] += q[obj_idx]*r;
            }
            if (g_flag!=Flag::POTENTIAL) {
                double r3 = r*r*r;
                ex[obj_idx] += q[src_idx]*dx*r3;
                ey[obj_idx] += q[src_idx]*dy*r3;
                ez[obj_idx] += q[src_idx]*dz*r3;
                ex[src_idx] -= q[obj_idx]*dx*r3;
                ey[src_idx] -= q[obj_idx]*dy*r3;
                ez[src_idx] -= q[obj_idx]*dz*r3;
            }
            src_idx = g_ptclist[src_idx];
        }
        obj_idx = g_ptclist[obj_idx];
    }
    return 0;
}

////Convert the multipole expansion into local expansion inside each other for two well separated boxes
int well_separated(unsigned long int obj_idx, Box &obj_box, unsigned long int src_idx, Box &src_box, double boxsize,
                   double *nabla_r) {
    //use the first n elements in g_local_expns and g_multipole_expns as g_scratch
    multipole_to_local(&g_multipole_expns[obj_idx*TOTAL_ELEMENT_NUMBER], obj_box.center[0], obj_box.center[1],
                       obj_box.center[2],src_box.center[0],src_box.center[1],src_box.center[2], nabla_r, boxsize,
                       g_local_expns);
    multipole_to_local(&g_multipole_expns[src_idx*TOTAL_ELEMENT_NUMBER], src_box.center[0], src_box.center[1],
                       src_box.center[2],obj_box.center[0],obj_box.center[1],obj_box.center[2], nabla_r, boxsize,
                       g_multipole_expns);
    for (int i=0; i<TOTAL_ELEMENT_NUMBER; ++i) {
        g_local_expns[src_idx*TOTAL_ELEMENT_NUMBER+i] += g_local_expns[i];
        g_local_expns[obj_idx*TOTAL_ELEMENT_NUMBER+i] += g_multipole_expns[i];
    }
    return 0;
}

////Translate the local expansion from the parent box to its child box
int local_exp_from_parent(unsigned long int parent_idx, Box &parent_box, unsigned long int child_idx, Box &child_box,
                          double *rho_tensor) {
    //use the first n elements in g_local_expns as g_scratch
    local_to_local(&g_local_expns[parent_idx*TOTAL_ELEMENT_NUMBER], parent_box.center[0], parent_box.center[1],
                   parent_box.center[2], child_box.center[0], child_box.center[1], child_box.center[2], rho_tensor,
                   g_local_expns);
    for(int i=0; i<TOTAL_ELEMENT_NUMBER; ++i) {
        g_local_expns[child_idx*TOTAL_ELEMENT_NUMBER+i] += g_local_expns[i];
    }
    return 0;
}
//
////Calculate the potential inside a childless box
int childless_box_potential(unsigned long int box_idx, Box &box, double *x, double *y, double *z, double *q,
                            double *phi, double *ex, double *ey, double *ez) {
    unsigned long int ptc_idx = box.first_ptcl;
    while (ptc_idx<TOTAL_PARTICLE_NUMBER) {
        double null = 0;
        switch (g_flag) {
            case Flag::POTENTIAL: {
                local_to_charge(&g_local_expns[box_idx*TOTAL_ELEMENT_NUMBER], box.center[0], box.center[1],
                                box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx], phi[ptc_idx], null, null, null);
                break;
            }
            case Flag::FIELD: {
                local_to_charge(&g_local_expns[box_idx*TOTAL_ELEMENT_NUMBER], box.center[0], box.center[1],
                                box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx], null,
                                ex[ptc_idx], ey[ptc_idx], ez[ptc_idx]);
                break;
            }
            case Flag::BOTH: {
                local_to_charge(&g_local_expns[box_idx*TOTAL_ELEMENT_NUMBER], box.center[0], box.center[1],
                                box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx],
                                phi[ptc_idx], ex[ptc_idx], ey[ptc_idx], ez[ptc_idx]);
                break;
            }
            default: {
                assert(false);
            }
        }
        ptc_idx = g_ptclist[ptc_idx];
    }
    return 0;
}


//Two boxes are ill separated and the larger one is childless
int ill_separated(unsigned long int large_idx, Box &large_box, unsigned long int small_idx, Box &small_box,
                  double *x, double *y, double *z, double *q, double *phi, double *ex, double *ey, double *ez) {
    //Calculate the potential inside the large box from the multipole expansion of the small box
    unsigned long int ptc_idx = large_box.first_ptcl;
    while (ptc_idx<TOTAL_PARTICLE_NUMBER) {
        double null = 0;
        switch (g_flag) {
            case Flag::POTENTIAL: {
                multipole_to_charge(&g_multipole_expns[small_idx*TOTAL_ELEMENT_NUMBER], small_box.center[0],
                                    small_box.center[1], small_box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx],
                                    phi[ptc_idx], null, null, null);
                break;
            }
            case Flag::FIELD: {
                multipole_to_charge(&g_multipole_expns[small_idx*TOTAL_ELEMENT_NUMBER], small_box.center[0],
                                    small_box.center[1], small_box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx],
                                    null, ex[ptc_idx], ey[ptc_idx], ez[ptc_idx]);
                break;
            }
            case Flag::BOTH: {
                multipole_to_charge(&g_multipole_expns[small_idx*TOTAL_ELEMENT_NUMBER], small_box.center[0],
                                    small_box.center[1], small_box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx],
                                    phi[ptc_idx], ex[ptc_idx], ey[ptc_idx], ez[ptc_idx]);
                break;
            }
            default: {
                assert(false);
            }
        }
        ptc_idx = g_ptclist[ptc_idx];
    }
    //Calculate the local expansion inside the small box using the charges inside the large box
    charge_to_local_traceless(large_box,q,x,y,z,small_box.center[0],small_box.center[1],small_box.center[2],
                              g_local_expns);
    for(int i=0; i<TOTAL_ELEMENT_NUMBER;++i){
        g_local_expns[small_idx*TOTAL_ELEMENT_NUMBER+i] += g_local_expns[i];
    }
    return 0;
}

//Check the descent of the parent box
int check_descent(vector<Box> &tree, unsigned long int childless_idx, Box &childless_box, Box &parent_box, double *x,
                  double *y, double *z, double *q, double *phi, double *ex, double *ey, double *ez) {
    //Create a stack of the boxes to check
    vector<unsigned long int> box_stack;
    for(int i=0;i<parent_box.n_child;++i) box_stack.push_back(parent_box.child[i]);

    unsigned long int check_idx;    //index of the box to check
    while (!box_stack.empty()) {
        check_idx = box_stack.back();
        box_stack.pop_back();
        Box & check_box = tree[check_idx];
        switch (separate(childless_box,check_box)) {
            case 0:{    //adjacent
                if (check_box.n_child==0) { //Check box is also childless, use Coulomb formula to calculate the potential
                    coulomb_potential(childless_box,check_box, x, y, z, q, phi, ex, ey, ez);
                }
                else{   //Check box if parent, add its children into stack
                    for (int i=0; i<check_box.n_child; ++i)  box_stack.push_back(check_box.child[i]);
                }
                break;
            }
            case 1:{    //ill separated
                ill_separated(childless_idx, childless_box, check_idx, check_box, x, y, z, q, phi, ex, ey, ez);
                break;
            }
            default:{
                cout<<"Warning, unexpected relation between the boxes in check_descent()"<<endl;
            }
        }
    }
    return 0;
}

//check the relation between box itr and the descent of box clg_idx
int check_colleague_child(vector<Box> &tree, unsigned long int itr, unsigned long int clg_idx, double &current_boxsize,
                          double *nabla_r, double *x, double *y, double *z, double *q, double *phi,
                          double *ex, double *ey, double *ez) {
    for(int j=0; j<tree[clg_idx].n_child;++j){
        unsigned long int clg_child_idx = tree[clg_idx].child[j];
        if (clg_child_idx>itr) {
            switch (separate(tree[itr],tree[clg_child_idx])) {
            case 0:{    //These two boxes are adjacent.
                if (tree[itr].n_child==0) {   //box itr is childless
                        if (tree[clg_child_idx].n_child==0) {     //box clg_child_idx is childless
                            coulomb_potential(tree[itr],tree[clg_child_idx], x, y, z, q, phi, ex, ey, ez);
                        }
                        else {       //child clg_child_idx is NOT childless
                            //check the descent of box clg_child_idx
                            check_descent(tree, itr, tree[itr], tree[clg_child_idx], x, y, z, q, phi, ex, ey, ez);
                        }
                }
                else {
                        if(tree[clg_child_idx].n_child==0){     //box clg_child_idx is childless
                            //check the descent of box itr
                            check_descent(tree, clg_child_idx, tree[clg_child_idx], tree[itr], x, y, z, q,
                                          phi, ex, ey, ez);
                        }
                }
                break;
            }
            case 2:{    //These two boxes are well separated.
                if (tree[itr].box_size<0.75*current_boxsize) {
                    update_nabla_r(nabla_r);
                    current_boxsize = tree[itr].box_size;
                }
                well_separated(itr, tree[itr], clg_child_idx, tree[clg_child_idx], current_boxsize, nabla_r);
                break;
            }
            default:
                cout<<"Warning, unexpected relation between the boxes in check_colleague_child()"<<endl;
            }
        }
    }
    return 0;
}

//Calculate the Coulomb potential on each particle by MLFMA
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, unsigned int max_rank, unsigned int n_ptc_box,
        double *phi, double *ex, double *ey, double *ez) {

    while (n_ptc_box+1>n_ptc) n_ptc_box /= 2;
    if (n_ptc_box<1) n_ptc_box = 1;
    vector<Box> tree;
	vector<Colleague> clg;
    //initialize g_ptclist, which records the particles in each childless box
	g_ptclist = new unsigned long int[n_ptc];
	memset(g_ptclist, 0, n_ptc*sizeof(unsigned long int));
	create_tree(x,y,z,n_ptc, n_ptc_box, tree, g_ptclist);
    configure_fmm(max_rank, n_ptc, tree.size());
	create_colleague(tree, clg);
    //Save the coefs for nabla operator
    g_nabla_coef = new double[kRankNNablaLength[MAX_RANK]];
    calc_nabla_1_emement_coef(g_nabla_coef);

    calc_multipole(tree, q, x, y, z);

    //set zero for output potential
    if(g_flag!=Flag::FIELD) memset(phi, 0, n_ptc*sizeof(double));
    if(g_flag!=Flag::POTENTIAL) {
        memset(ex, 0, n_ptc*sizeof(double));
        memset(ey, 0, n_ptc*sizeof(double));
        memset(ez, 0, n_ptc*sizeof(double));
    }
    double * rho_tensor = new double[8*TOTAL_ELEMENT_NUMBER];
    double current_parentbox_size = 0.5*tree[0].box_size;
    calc_rho_tensor(current_parentbox_size,rho_tensor);

    double * nabla_r = new double[16*TOTAL_ELEMENT_NUMBER];
    double  current_boxsize = 0.25*tree[0].box_size;
    calc_nabla_r(current_boxsize,nabla_r);

    for (unsigned long int itr=1; itr<tree.size();++itr) {
        unsigned long int parent_idx = tree[itr].parent;
        if (parent_idx>0) {
            if (tree[parent_idx].box_size<0.75*current_parentbox_size) {
                    update_rho_tensor(rho_tensor);
                    current_parentbox_size = tree[parent_idx].box_size;
            }
            //Inherit the local expansion from the parent box
            local_exp_from_parent(parent_idx, tree[parent_idx], itr, tree[itr], rho_tensor);
            //Check all the child boxes of the colleagues of b's parent box
            for (int i=0; clg[parent_idx].clg[i]>0;++i) {
                unsigned long int clg_idx = clg[parent_idx].clg[i];
                check_colleague_child(tree, itr, clg_idx, current_boxsize, nabla_r, x, y, z, q, phi, ex, ey, ez);
            }
        }
        else {
            check_colleague_child(tree, itr, parent_idx, current_boxsize, nabla_r, x, y, z, q, phi, ex, ey, ez);
        }
        if (tree[itr].n_child==0) {
            coulomb_potential(tree[itr], x, y, z, q, phi, ex, ey, ez);
            childless_box_potential(itr, tree[itr], x, y, z, q, phi, ex, ey, ez);
        }
    }
    delete[] g_ptclist;
    delete[] rho_tensor;
    delete[] nabla_r;
    delete[] g_nabla_coef;
	end_fmm();
	return 0;
}

//Calculate the Coulomb field on each particle by MLFMA
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, unsigned int max_rank, unsigned int n_ptc_box,
        double *ex, double *ey, double *ez) {
    g_flag = Flag::FIELD;
    fmm(x, y, z, q, n_ptc, max_rank, n_ptc_box, nullptr, ex, ey, ez);
    g_flag = Flag::BOTH;
    return 0;
}

//Calculate the Coulomb potential on each particle by MLFMA
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, unsigned int max_rank, unsigned int n_ptc_box,
        double *phi) {
    g_flag = Flag::POTENTIAL;
    fmm(x, y, z, q, n_ptc, max_rank, n_ptc_box, phi, nullptr, nullptr, nullptr);
    g_flag = Flag::BOTH;
    return 0;
}
