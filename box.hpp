/**********************************
box.hpp
Define the box data structure and the colleague list for each box
Declare functions to manipulate the box structure

version 1.0
By He Zhang, 02/2015

***********************************/


#ifndef BOX_HPP
#define BOX_HPP

#include <vector>
#include <cstdio>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <cmath>

using std::vector;
using std::cout;
using std::endl;

typedef struct Box{
	double center[3] = {0,0,0};
	unsigned long int parent = 0;
	unsigned long int child[8] = {0};
	unsigned long int first_ptcl = 0;
	unsigned int n_ptcl = 0;
	int n_child = 0;
	double box_size = 0;
} Box;

typedef struct Colleague{
	unsigned long int clg[28] = {0}; //At most 27 colleagues, the last 0 means the end.
} Colleague;

//Create the hierarchical tree structure of the boxes
int create_tree(double * x, double * y, double * z, const unsigned long int n, const unsigned int s, vector<Box> &tree, unsigned long int * list);
//Create the colleague list for all boxes
int create_colleague(vector<Box> &box, vector<Colleague> &clg);
//Find the relation between two boxes: adjcent (return 0), ill separated (return 1), or well separated (return 2)
int separate(Box &box1, Box &box2);

std::ostream& operator<<(std::ostream& os, Box& box);
std::ostream& operator<<(std::ostream& os, Colleague& clg);
#endif
