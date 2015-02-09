#include "head.hpp"

//Calculate the multipole expansions for all boxes except for the root box
int calc_multipole(vector<Box> &tree, double * q, double * x, double * y, double * z){
    //use global ptclist, multipole_expns.

    for(unsigned long int itr=tree.size()-1; itr>0; --itr){
		if (tree[itr].n_child>0){  //parent box, take summation of the multipole expansion of the child boxes
                for (int i=0; i<tree[itr].n_child; ++i){
                        unsigned long int child_ptr = tree[itr].child[i];
                        //Translate the mltp of a child box to the parent box, use the first Number_of_total_element number of multipoles as scratch variable
                        Multipole_to_Multipole(tree[child_ptr].center[0],tree[child_ptr].center[1],tree[child_ptr].center[2],tree[itr].center[0],tree[itr].center[1],tree[itr].center[2], &multipole_expns[child_ptr*Number_of_total_element], multipole_expns );
                        for (int j=0; j<Number_of_total_element; ++j){
                            multipole_expns[itr*Number_of_total_element+j] += multipole_expns[j];
                        }
                }
		}
		else{	//childless box, calculate the multipole expansion from the charges
//            Box & box = tree[itr];
//            Charge_to_Multipole(box, q, x, y, z, ptclist, &multipole_expns[itr*Number_of_total_element]);
            Charge_to_Multipole(tree[itr], q, x, y, z, ptclist, &multipole_expns[itr*Number_of_total_element]);
		}
	}

	return 0;
}

//Calculate the potential using Coulomb formula in a childless box
int Coulomb_potential(Box &box, double * x, double * y, double * z, double * q, double * phi){
    //use global ptclist.

    unsigned long int obj_idx, src_idx;
    double r = 0; //store 1/r
    obj_idx = box.first_ptcl;
    while (obj_idx<Number_of_particle){
        src_idx = ptclist[obj_idx];
        while (src_idx<Number_of_particle){
            r = 1.0/sqrt(pow(x[obj_idx]-x[src_idx],2)+pow(y[obj_idx]-y[src_idx],2)+pow(z[obj_idx]-z[src_idx],2));
            phi[obj_idx] += q[src_idx]*r;
            phi[src_idx] += q[obj_idx]*r;
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
    double r = 0; //store 1/r
    obj_idx = obj_box.first_ptcl;
    while (obj_idx<Number_of_particle){
        src_idx = src_box.first_ptcl;
        while (src_idx<Number_of_particle){
            r = 1.0/sqrt(pow(x[obj_idx]-x[src_idx],2)+pow(y[obj_idx]-y[src_idx],2)+pow(z[obj_idx]-z[src_idx],2));
            phi[obj_idx] += q[src_idx]*r;
            phi[src_idx] += q[obj_idx]*r;
            src_idx = ptclist[src_idx];
        }
        obj_idx = ptclist[obj_idx];
    }
    return 0;
}

//Convert the multipole expansion into local expansion inside each other for two well separated boxes
int well_seperated(unsigned long int obj_idx, Box &obj_box, unsigned long int src_idx, Box &src_box){
    //use global local_expns, multipole_expns.

    //use the first n elements in local_expns and multipole_expns as scratch
    Multipole_to_Local(&multipole_expns[obj_idx*Number_of_total_element], obj_box.center[0], obj_box.center[1], obj_box.center[2],src_box.center[0],src_box.center[1],src_box.center[2],local_expns);
    Multipole_to_Local(&multipole_expns[src_idx*Number_of_total_element], src_box.center[0], src_box.center[1], src_box.center[2],obj_box.center[0],obj_box.center[1],obj_box.center[2],multipole_expns);

    for (int i=0; i<Number_of_total_element; ++i){
        local_expns[obj_idx*Number_of_total_element+i] += local_expns[i];
        local_expns[src_idx*Number_of_total_element+i] += multipole_expns[i];
    }
    return 0;
}

//Transfer the local expansion to the child box from its parent box
int local_exp_from_parent(unsigned long int parent_idx, Box &parent_box, unsigned long int child_idx, Box &child_box){
    //use global local_expns.

    //use the first n elements in local_expns as scratch
    Local_to_Local(&local_expns[parent_idx*Number_of_total_element], parent_box.center[0], parent_box.center[1], parent_box.center[2], child_box.center[0], child_box.center[1], child_box.center[2], local_expns);
    for(int i=0; i<Number_of_total_element; ++i){
        local_expns[child_idx*Number_of_total_element+i] += local_expns[i];
    }
    return 0;
}

//Calculate the potential inside a childless box
int Childless_box_potential(unsigned long int box_idx, Box &box, double * x, double * y, double * z, double * q, double * phi){
    //use global local_expns, ptclist.

    unsigned long int ptc_idx = box.first_ptcl;
    while(ptc_idx<Number_of_particle){
        phi[ptc_idx] += LocalPotential(&local_expns[box_idx*Number_of_total_element], box.center[0], box.center[1], box.center[2], x[ptc_idx], y[ptc_idx], z[ptc_idx]);
        ptc_idx = ptclist[ptc_idx];
    }
    return 0;
}

int fmm(double * x, double * y, double * z, double * q, unsigned long int n_ptc, int max_rank, int n_ptc_box, double * phi){

    vector<Box> tree;
	vector<Colleague> clg;

    //initialize ptclist, which records the particles in each childless box
	ptclist = new unsigned long int[n_ptc];
	memset(ptclist, 0, n_ptc*sizeof(unsigned long int));

	create_tree(x,y,z,n_ptc, n_ptc_box, tree, ptclist);

    configure_fmm(max_rank, n_ptc, tree.size());
	cout<<n_Max_rank<<' '<<Number_of_total_element<<endl;


//	for(unsigned long int i=0; i<n_ptc; ++i){
//		cout<<i<<' '<<ptclist[i]<<' '<<x[i]<<' '<<y[i]<<' '<<z[i]<<endl;
//	}
//
//	int i=0;
//	for(auto itr=tree.begin(); itr!=tree.end(); ++itr){
//		cout<<i<<' '<<*itr<<endl;
//		++i;
//	}
//
//	cout<<clg.size()<<endl;
//
	create_colleague(tree, clg);
//	cout<<clg.size()<<endl;
//
//
//	cout<<endl;
//	for(unsigned long int i=0; i<clg.size(); ++i){
//		cout<<clg[i];
//	}

    calc_multipole(tree, q, x, y, z);

    //set zero for output potential
    memset(phi, 0, n_ptc*sizeof(double));

    for(unsigned long int itr=1; itr<tree.size();++itr){
        unsigned long int parent_idx = tree[itr].parent;
        if(tree[itr].parent>0)  local_exp_from_parent(parent_idx, tree[parent_idx], itr, tree[itr] );   //Inherit the local expansion from the parent box
        for(int i=0; clg[parent_idx].clg[i]>0;++i){     //Check all the child boxes of the colleagues of b's parent box
            unsigned long int clg_idx = clg[parent_idx].clg[i];
            for(int j=0; j<tree[clg_idx].n_child;++j){
                unsigned long int clg_child_idx = tree[clg_idx].child[j];
                if(clg_child_idx>itr){
                    switch(separate(tree[itr],tree[clg_child_idx])){
                    case 0:{    //These two boxes are adjacent.
                        if(tree[itr].n_child==0){   //box itr is childless
                                Coulomb_potential(tree[itr],x,y,z,q,phi);
                                if(tree[clg_child_idx].n_child==0){     //box clg_child_idx is childless
                                    Coulomb_potential(tree[itr],tree[clg_child_idx],x,y,z,q,phi);
                                }
                                else{       //child clg_child_idx is NOT childless
                                        //check the descent of box clg_child_idx

                                }
                        }
                        else{
                                if(tree[clg_child_idx].n_child==0){     //box clg_child_idx is childless
                                    //check the descent of box itr
                                }
                        }
                        break;
                    }
                    case 2:{    //These two boxes are well separated.
                        well_seperated(itr, tree[itr], clg_child_idx, tree[clg_child_idx]);
                        break;
                    }
                    default:
                        cout<<"Warning, unexpected relation between the boxes"<<endl;
                    }
                }
            }
        }

        if(tree[itr].n_child==0) Childless_box_potential(itr, tree[itr], x, y, z, q, phi);

    }

//local_exp_from_parent(unsigned long int parent_idx, Box &parent_box, unsigned long int child_idx, Box &child_box)
	end_fmm();


	return 0;
}
