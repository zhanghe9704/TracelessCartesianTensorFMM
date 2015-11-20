/**********************************
FMMframe.cpp

Frame of the MLFMA
Define the function to calculate the Coulomb potential using MLFMA and all related functions used in the frame

version 2.0
By He Zhang, 04/06/2015
Calculate 3D field

version 1.0
By He Zhang, 02/2015

***********************************/

#include <cassert>
#include "global.h"

//Calculate the multipole expansions for all boxes except for the root box
int calc_multipole(vector<Box> &tree, double * q, double * x, double * y, double * z){
    //use global ptclist, multipole_expns.

    double * multipole_coef = new double[multipole_coef_length[n_Max_rank]];
    double current_childbox_size = tree[tree.size()-1].box_size;
    multipole_to_multipole_coef(current_childbox_size, multipole_coef);

    for(unsigned long int itr=tree.size()-1; itr>0; --itr){
		if (tree[itr].n_child>0){  //parent box, take summation of the multipole expansion of the child boxes
                for (int i=0; i<tree[itr].n_child; ++i){
                        unsigned long int child_ptr = tree[itr].child[i];
////                        Translate the mltp of a child box to the parent box, use the first Number_of_total_element number of multipoles as scratch variable
                        if (tree[child_ptr].box_size>1.5*current_childbox_size){
                            update_multipole_to_multipole_coef(multipole_coef);
                            current_childbox_size = tree[child_ptr].box_size;
                        }
                        Multipole_to_Multipole(tree[child_ptr].center[0],tree[child_ptr].center[1],tree[child_ptr].center[2],tree[itr].center[0],tree[itr].center[1],tree[itr].center[2], multipole_coef, &multipole_expns[child_ptr*Number_of_total_element], multipole_expns );
                        for (int j=0; j<Number_of_total_element; ++j){
                            multipole_expns[itr*Number_of_total_element+j] += multipole_expns[j];
                        }
                }
		}
		else{	//childless box, calculate the multipole expansion from the charges
            Charge_to_Multipole(tree[itr], q, x, y, z, ptclist, &multipole_expns[itr*Number_of_total_element]);
		}
	}
    delete[] multipole_coef;
	return 0;
}

//Calculate the potential using Coulomb formula in a childless box
int Coulomb_potential(Box &box, double * x, double * y, double * z, double * q, double * phi){
    //use global ptclist.

    unsigned long int obj_idx, src_idx;
    obj_idx = box.first_ptcl;
    while (obj_idx<Number_of_particle){
        src_idx = ptclist[obj_idx];
        while (src_idx<Number_of_particle){
            double dx = x[obj_idx]-x[src_idx];
            double dy = y[obj_idx]-y[src_idx];
            double dz = z[obj_idx]-z[src_idx];
            double r = 1/sqrt(dx*dx+dy*dy+dz*dz);
            phi[obj_idx] += q[src_idx]*r;
            phi[src_idx] += q[obj_idx]*r;
            src_idx = ptclist[src_idx];
        }
        obj_idx = ptclist[obj_idx];
    }
    return 0;
}

int coulomb_potential(Box &box, double *x, double *y, double *z, double *q, double *phi,
                      double *ex, double *ey, double *ez){
    //use global ptclist.

    unsigned long int obj_idx, src_idx;
    obj_idx = box.first_ptcl;
    while (obj_idx<Number_of_particle){
        src_idx = ptclist[obj_idx];
        while (src_idx<Number_of_particle){
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
            src_idx = ptclist[src_idx];
        }
        obj_idx = ptclist[obj_idx];
    }
    return 0;
}

//Calculate the potential using Coulomb formula in two adjacent childless box
int Coulomb_potential(Box &obj_box, Box &src_box, double * x, double * y, double * z, double * q, double * phi){
    //use global ptclist.

    unsigned long int obj_idx, src_idx;
    obj_idx = obj_box.first_ptcl;
    while (obj_idx<Number_of_particle){
        src_idx = src_box.first_ptcl;
        while (src_idx<Number_of_particle){
            double dx = x[obj_idx]-x[src_idx];
            double dy = y[obj_idx]-y[src_idx];
            double dz = z[obj_idx]-z[src_idx];
            double r = 1/sqrt(dx*dx+dy*dy+dz*dz);
            phi[obj_idx] += q[src_idx]*r;
            phi[src_idx] += q[obj_idx]*r;
            src_idx = ptclist[src_idx];
        }
        obj_idx = ptclist[obj_idx];
    }
    return 0;
}

int coulomb_potential(Box &obj_box, Box &src_box, double *x, double *y, double *z, double *q,
                      double *phi, double *ex, double *ey, double *ez){
    //use global ptclist.

    unsigned long int obj_idx, src_idx;
    obj_idx = obj_box.first_ptcl;
    while (obj_idx<Number_of_particle){
        src_idx = src_box.first_ptcl;
        while (src_idx<Number_of_particle){
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
            src_idx = ptclist[src_idx];
        }
        obj_idx = ptclist[obj_idx];
    }
    return 0;
}
//
////Convert the multipole expansion into local expansion inside each other for two well separated boxes
int well_separated(unsigned long int obj_idx, Box &obj_box, unsigned long int src_idx, Box &src_box, double boxsize, double * Nabla_R){
    //use global local_expns, multipole_expns.

    //use the first n elements in local_expns and multipole_expns as scratch
    Multipole_to_Local(&multipole_expns[obj_idx*Number_of_total_element], obj_box.center[0], obj_box.center[1], obj_box.center[2],src_box.center[0],src_box.center[1],src_box.center[2], Nabla_R, boxsize, local_expns);
    Multipole_to_Local(&multipole_expns[src_idx*Number_of_total_element], src_box.center[0], src_box.center[1], src_box.center[2],obj_box.center[0],obj_box.center[1],obj_box.center[2], Nabla_R, boxsize, multipole_expns);

    for (int i=0; i<Number_of_total_element; ++i){
        local_expns[src_idx*Number_of_total_element+i] += local_expns[i];
        local_expns[obj_idx*Number_of_total_element+i] += multipole_expns[i];
    }
    return 0;
}

////Translate the local expansion from the parent box to its child box
int local_exp_from_parent(unsigned long int parent_idx, Box &parent_box, unsigned long int child_idx, Box &child_box, double * Rho_Tensor){
    //use global local_expns.

    //use the first n elements in local_expns as scratch
    Local_to_Local(&local_expns[parent_idx*Number_of_total_element], parent_box.center[0], parent_box.center[1], parent_box.center[2], child_box.center[0], child_box.center[1], child_box.center[2], Rho_Tensor, local_expns);
//    Local_to_Local(&local_expns[parent_idx*Number_of_total_element], parent_box.center[0], parent_box.center[1], parent_box.center[2], child_box.center[0], child_box.center[1], child_box.center[2], local_expns);
    for(int i=0; i<Number_of_total_element; ++i){
        local_expns[child_idx*Number_of_total_element+i] += local_expns[i];
    }
    return 0;
}
//
////Calculate the potential inside a childless box
int Childless_box_potential(unsigned long int box_idx, Box &box, double *x, double *y, double *z, double *q,
                            double *phi, double *ex, double *ey, double *ez){
    //use global local_expns, ptclist.

    unsigned long int ptc_idx = box.first_ptcl;
    while(ptc_idx<Number_of_particle){
        double null = 0;
//        phi[ptc_idx] += LocalPotential(&local_expns[box_idx*Number_of_total_element], box.center[0], box.center[1], box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx]);
        switch (g_flag){
            case Flag::POTENTIAL: {
                local_to_charge(&local_expns[box_idx*Number_of_total_element], box.center[0], box.center[1],
                                box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx], phi[ptc_idx], null, null, null);
                break;
            }
            case Flag::FIELD: {
                local_to_charge(&local_expns[box_idx*Number_of_total_element], box.center[0], box.center[1],
                                box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx], null,
                                ex[ptc_idx], ey[ptc_idx], ez[ptc_idx]);
                break;
            }
            case Flag::BOTH: {
                local_to_charge(&local_expns[box_idx*Number_of_total_element], box.center[0], box.center[1],
                                box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx],
                                phi[ptc_idx], ex[ptc_idx], ey[ptc_idx], ez[ptc_idx]);
                break;
            }
            default: {
                assert(false);
            }
        }
        ptc_idx = ptclist[ptc_idx];
    }
    return 0;
}


//Two boxes are ill separated and the larger one is childless
int ill_separated(unsigned long int large_idx, Box &large_box, unsigned long int small_idx, Box &small_box,
                  double *x, double *y, double *z, double *q, double *phi, double *ex, double *ey, double *ez){
    //use global multipole_expns, local_expns, Number_of_particle, Number_of_total_element, ptclist.

    //Calculate the potential inside the large box from the multipole expansion of the small box
    unsigned long int ptc_idx = large_box.first_ptcl;
    while(ptc_idx<Number_of_particle){
//        phi[ptc_idx] += MultipolePotential(&multipole_expns[small_idx*Number_of_total_element], small_box.center[0], small_box.center[1], small_box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx]);
        double null = 0;
        switch(g_flag){
            case Flag::POTENTIAL: {
                multipole_to_charge(&multipole_expns[small_idx*Number_of_total_element], small_box.center[0],
                                    small_box.center[1], small_box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx],
                                    phi[ptc_idx], null, null, null);
                break;
            }
            case Flag::FIELD: {
                multipole_to_charge(&multipole_expns[small_idx*Number_of_total_element], small_box.center[0],
                                    small_box.center[1], small_box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx],
                                    null, ex[ptc_idx], ey[ptc_idx], ez[ptc_idx]);
                break;
            }
            case Flag::BOTH: {
                multipole_to_charge(&multipole_expns[small_idx*Number_of_total_element], small_box.center[0],
                                    small_box.center[1], small_box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx],
                                    phi[ptc_idx], ex[ptc_idx], ey[ptc_idx], ez[ptc_idx]);
                break;
            }
            default: {
                assert(false);
            }
        }
        ptc_idx = ptclist[ptc_idx];
    }

    //Calculate the local expansion inside the small box using the charges inside the large box
    Charge_to_Local_traceless(large_box,ptclist,q,x,y,z,small_box.center[0],small_box.center[1],small_box.center[2],local_expns);
    for(int i=0; i<Number_of_total_element;++i){
        local_expns[small_idx*Number_of_total_element+i] += local_expns[i];
    }

    return 0;
}

//Check the descent of the parent box
int check_descent(vector<Box> &tree, unsigned long int childless_idx, Box &childless_box, Box &parent_box, double *x, double *y, double *z, double *q,
                  double *phi, double *ex, double *ey, double *ez){

    //Create a stack of the boxes to check
    vector<unsigned long int> box_stack;
    for(int i=0;i<parent_box.n_child;++i) box_stack.push_back(parent_box.child[i]);

    unsigned long int check_idx;    //index of the box to check
    while(!box_stack.empty()){
        check_idx = box_stack.back();
        box_stack.pop_back();
        Box & check_box = tree[check_idx];
        switch(separate(childless_box,check_box)){
            case 0:{    //adjacent
                if(check_box.n_child==0){   //Check box is also childless, use Coulomb formula to calculate the potential
                    coulomb_potential(childless_box,check_box, x, y, z, q, phi, ex, ey, ez);
                }
                else{   //Check box if parent, add its children into stack
                    for(int i=0; i<check_box.n_child; ++i)  box_stack.push_back(check_box.child[i]);
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
                          double *Nabla_R, double *x, double *y, double *z, double *q, double *phi,
                          double *ex, double *ey, double *ez){

            for(int j=0; j<tree[clg_idx].n_child;++j){
                unsigned long int clg_child_idx = tree[clg_idx].child[j];
                if(clg_child_idx>itr){
                    switch(separate(tree[itr],tree[clg_child_idx])){
                    case 0:{    //These two boxes are adjacent.
                        if(tree[itr].n_child==0){   //box itr is childless
                                //Coulomb_potential(tree[itr],x,y,z,q,phi);
                                if(tree[clg_child_idx].n_child==0){     //box clg_child_idx is childless
                                    coulomb_potential(tree[itr],tree[clg_child_idx], x, y, z, q, phi, ex, ey, ez);
                                }
                                else{       //child clg_child_idx is NOT childless
                                    //check the descent of box clg_child_idx
                                    check_descent(tree, itr, tree[itr], tree[clg_child_idx], x, y, z, q, phi, ex, ey, ez);
                                }
                        }
                        else{
                                if(tree[clg_child_idx].n_child==0){     //box clg_child_idx is childless
                                    //check the descent of box itr
                                    check_descent(tree, clg_child_idx, tree[clg_child_idx], tree[itr], x, y, z, q,
                                                  phi, ex, ey, ez);
                                }
                        }
                        break;
                    }
                    case 2:{    //These two boxes are well separated.
                        if(tree[itr].box_size<0.75*current_boxsize){
                            update_Nabla_R(Nabla_R);
                            current_boxsize = tree[itr].box_size;
                        }
//                        well_separated(itr, tree[itr], clg_child_idx, tree[clg_child_idx]);
                        well_separated(itr, tree[itr], clg_child_idx, tree[clg_child_idx], current_boxsize, Nabla_R);
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
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double *phi, double *ex, double *ey, double *ez){

    vector<Box> tree;
	vector<Colleague> clg;

    //initialize ptclist, which records the particles in each childless box
	ptclist = new unsigned long int[n_ptc];
	memset(ptclist, 0, n_ptc*sizeof(unsigned long int));

	create_tree(x,y,z,n_ptc, n_ptc_box, tree, ptclist);
    configure_fmm(max_rank, n_ptc, tree.size());
	create_colleague(tree, clg);

    //Save the coefs for Nabla operator
    Nabla_1_element_r_coef = new double[Nabla_1_element_r_length[n_Max_rank]];
    Calc_Nabla_1_emement_coef(Nabla_1_element_r_coef);

    calc_multipole(tree, q, x, y, z);

    //set zero for output potential
    if(g_flag!=Flag::FIELD) memset(phi, 0, n_ptc*sizeof(double));
    if(g_flag!=Flag::POTENTIAL) {
        memset(ex, 0, n_ptc*sizeof(double));
        memset(ey, 0, n_ptc*sizeof(double));
        memset(ez, 0, n_ptc*sizeof(double));
    }

    double * Rho_Tensor = new double[8*Number_of_total_element];
    double current_parentbox_size = 0.5*tree[0].box_size;
    Calc_Rho_Tensor(current_parentbox_size,Rho_Tensor);

    double * Nabla_R = new double[16*Number_of_total_element];
    double  current_boxsize = 0.25*tree[0].box_size;
    Calc_Nabla_R(current_boxsize,Nabla_R);


    for(unsigned long int itr=1; itr<tree.size();++itr){
        unsigned long int parent_idx = tree[itr].parent;

        if(parent_idx>0){
            if(tree[parent_idx].box_size<0.75*current_parentbox_size) {
                    update_Rho_Tensor(Rho_Tensor);
                    current_parentbox_size = tree[parent_idx].box_size;
            }
            //Inherit the local expansion from the parent box
            local_exp_from_parent(parent_idx, tree[parent_idx], itr, tree[itr], Rho_Tensor);
            //Check all the child boxes of the colleagues of b's parent box
            for(int i=0; clg[parent_idx].clg[i]>0;++i){
                unsigned long int clg_idx = clg[parent_idx].clg[i];
                check_colleague_child(tree, itr, clg_idx, current_boxsize, Nabla_R, x, y, z, q, phi, ex, ey, ez);
            }
        }
        else{
            check_colleague_child(tree, itr, parent_idx, current_boxsize, Nabla_R, x, y, z, q, phi, ex, ey, ez);
        }

        if(tree[itr].n_child==0){
            coulomb_potential(tree[itr], x, y, z, q, phi, ex, ey, ez);
            Childless_box_potential(itr, tree[itr], x, y, z, q, phi, ex, ey, ez);
        }

    }

    delete[] ptclist;
    delete[] Rho_Tensor;
    delete[] Nabla_R;
    delete[] Nabla_1_element_r_coef;
	end_fmm();

	return 0;
}

//Calculate the Coulomb field on each particle by MLFMA
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, int max_rank, int n_ptc_box,
        double *ex, double *ey, double *ez) {
    g_flag = Flag::FIELD;
    return fmm(x, y, z, q, n_ptc, max_rank, n_ptc_box, nullptr, ex, ey, ez);
}

//Calculate the Coulomb potential on each particle by MLFMA
int fmm(double *x, double *y, double *z, double *q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double *phi) {
    g_flag = Flag::POTENTIAL;
    return fmm(x, y, z, q, n_ptc, max_rank, n_ptc_box, phi, nullptr, nullptr, nullptr);
}