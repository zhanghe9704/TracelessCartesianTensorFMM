#include "box.hpp"

int find_root_center(double *x, double * y, double * z, const unsigned long int N, double &cx, double &cy, double &cz, double &size){
	double max_x = x[0];
	double max_y = y[0];
	double max_z = z[0];
	double min_x = x[0];
	double min_y = y[0];
	double min_z = z[0];

	for(unsigned long int i=1; i<N; ++i){
		if (max_x<x[i]) max_x = x[i];
		if (max_y<y[i]) max_y = y[i];
		if (max_z<z[i]) max_z = z[i];
		if (min_x>x[i]) min_x = x[i];
		if (min_y>y[i]) min_y = y[i];
		if (min_z>z[i]) min_z = z[i];
	}

	cx = 0.5*(max_x+min_x);
	cy = 0.5*(max_y+min_y);
	cz = 0.5*(max_z+min_z);

	size = max_x-min_x;
	if (size<(max_y-min_y)) size = max_y - min_y;
	if (size<(max_z-min_z)) size = max_z - min_z;
	size *= (1+1e-6);

	return 0;
}

int create_tree(double * x, double * y, double * z, const unsigned long int n, const unsigned int s, vector<Box> &tree, unsigned long int * list){

	Box empty_box;

//	unsigned long int max_ptcl_in_box = n;
//	unsigned int n_box = 0;
	double cx, cy, cz, size;	//the center and the size of a box

	//create the root box
	tree.push_back(empty_box);
	find_root_center(x, y, z, n, cx, cy, cz, size);
	tree[0].center[0] = cx;
	tree[0].center[1] = cy;
	tree[0].center[2] = cz;
	tree[0].box_size = size;
	tree[0].n_ptcl = n;

	//prepare for the particle list
	for(unsigned long int i=0; i<n; ++i){
		list[i] = i+1;
	}

	unsigned long int tmp_idx[8] = {n, n, n, n, n, n, n, n};
	int nx, ny, nz;

	for(unsigned long int i=0; i!=tree.size();++i){
		if(tree[i].n_ptcl>s){
			//split the box

			//go through all the particles in the box
			unsigned long int idx=tree[i].first_ptcl;
			while(idx!=n){	//when idx==n, reached the end of the particle list
			//for (unsigned long int idx=tree[i].first_ptcl; idx<tree[i].n_ptcl; ++idx){
				nx = (x[idx]-tree[i].center[0]>0)?1:0;
				ny = (y[idx]-tree[i].center[1]>0)?1:0;
				nz = (z[idx]-tree[i].center[2]>0)?1:0;

				int idx_cld = 4*nx+2*ny+nz;

				if (tree[i].child[idx_cld]>0){	//child box has created.
					//update the child box
					unsigned long int pos_cld = tree[i].child[idx_cld];
					tree[pos_cld].n_ptcl += 1;

					//udpate the linked list of the particles
					list[tmp_idx[idx_cld]] = idx;
					tmp_idx[idx_cld] = idx;


				}
				else{	//child box has NOT created.

					//create the child box
					unsigned long int pos_cld = tree.size();
					tree.push_back(empty_box);
					tree[pos_cld].parent = i;
					tree[pos_cld].box_size = 0.5*tree[i].box_size;
					tree[pos_cld].first_ptcl = idx;
					tree[pos_cld].n_ptcl +=1;

					if (nx==0) nx=-1;
					if (ny==0) ny=-1;
					if (nz==0) nz=-1;
					tree[pos_cld].center[0] = tree[i].center[0]+0.5*nx*tree[pos_cld].box_size;
					tree[pos_cld].center[1] = tree[i].center[1]+0.5*ny*tree[pos_cld].box_size;
					tree[pos_cld].center[2] = tree[i].center[2]+0.5*nz*tree[pos_cld].box_size;

					tmp_idx[idx_cld] = idx;	//record the last particle in the box

					//update the parent box
					tree[i].child[idx_cld] = pos_cld;
					tree[i].n_child += 1;
					tree[i].first_ptcl = n;


				}
				idx = list[idx];
			}

			tree[i].n_ptcl = 0;
			for (int i=0; i<8; ++i) {
				if (tmp_idx[i]<n){
					list[tmp_idx[i]] = n;
					tmp_idx[i] = n;
				}

			}
		}

	}


	//Make sure the index of the child boxes stored in the first n_child elements of the child array.
	for(unsigned long int i=0; i!=tree.size(); ++i){
		if (tree[i].n_child>0){
			int k=7;
			for(int j=0; j<tree[i].n_child; ++j){
				if (tree[i].child[j]==0){
					for(;k>j;--k){
						if (tree[i].child[k]>0) {
							std::swap(tree[i].child[j], tree[i].child[k]);
							break;
						}
					}
				}
			}
		}

	}
	return 0;
}


//Check the two boxes are adjcent (return 0), ill separated (return 1), or well separated (return 2)
int separate(Box &box1, Box &box2){
	double large_size, small_size;
	double distance[3] = {0};

	large_size = box1.box_size;
	small_size = box2.box_size;

	if (large_size<small_size) std::swap(large_size, small_size);

	for(int i=0; i<3; ++i){
		distance[i] = abs(box1.center[i]-box2.center[i]);
	}

	large_size = 1.5*large_size;
	small_size = 1.5*small_size;
	if ((distance[0]>large_size)+(distance[1]>large_size)+(distance[2]>large_size)){
		return 2;
	}
	else if ((distance[0]>small_size)+(distance[1]>small_size)+(distance[2]>small_size)){
		return 1;
	}
	else{
		return 0;
	}

}

int create_colleague(vector<Box> &box, vector<Colleague> &clg){

	clg.reserve(box.size());
	Colleague empty_clg;
	for(unsigned long int i=0; i<box.size(); ++i) clg.push_back(empty_clg); //[i] = empty_clg;

	//Set the colleague list for all level one boxes to be the other level one boxes
	for (int i=0; i<box[0].n_child; ++i){
		int idx = box[0].child[i];
		memcpy(clg[idx].clg, box[0].child, box[0].n_child*sizeof(unsigned long int));
		std::swap(clg[idx].clg[0],clg[idx].clg[i]);

	}

	//For finer level boxes
	for(unsigned long int i=box[0].n_child+1; i<box.size(); ++i){

		clg[i].clg[0] = i; //The first colleague is the box itself
		unsigned long int prnt = box[i].parent;	//find parent box

		int count = 1;

		//Check the parent box
		for(int j=0; j<box[prnt].n_child; ++j){
			if (i!=box[prnt].child[j]) {
				clg[i].clg[count] = box[prnt].child[j];
				count += 1;
			}
		}

		//Check other colleagues of the parent box
		for(int j=1; clg[prnt].clg[j]!=0; ++j){
			unsigned long int prnt_clg = clg[prnt].clg[j];
			for (int k=0; k<box[prnt_clg].n_child; ++k){
				unsigned long int prnt_clg_cld = box[prnt_clg].child[k];
				if (separate(box[i], box[prnt_clg_cld])==0){
					clg[i].clg[count] = prnt_clg_cld;
					count += 1;
				}
			}

		}
	}

	return 0;
}


//output a box for debug
std::ostream& operator<<(std::ostream& os, Box& box){
	os<<"Box center: ";
	for(int i=0; i<3; ++i) os<<box.center[i]<<' ';
	os<<endl;

	os<<"Box size: " << box.box_size <<endl;

	os<<"Parent box: "<< box.parent <<endl;
	os<<"Number of child boxes: "<< box.n_child << endl;
	os<<"Child boxes: ";
	for(int i=0; i<box.n_child; ++i) {
//		if (box.child[i]>0) os<< box.child[i] <<' ';
		os<< box.child[i] <<' ';
	}
	os << endl;

	os<<"Particles: " << box.n_ptcl <<' ' << box.first_ptcl<<endl;
}


std::ostream& operator<<(std::ostream& os, Colleague& clg){
	os<<clg.clg[0]<<' ';
	for(int i=1; clg.clg[i]>0; ++i) os<<clg.clg[i]<<' ';
	os<<endl;
}

