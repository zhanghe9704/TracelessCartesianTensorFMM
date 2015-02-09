#ifndef BOX_HPP
#define BOX_HPP

#include <vector>
#include <cstdio>
#include <iostream>
#include <utility>
#include <cstdlib>

using std::vector;
using std::cout;
using std::endl;

typedef struct Box{
	double center[3] = {0,0,0};
	unsigned long int parent = 0;
	unsigned long int child[8] = {0};
	unsigned long int first_ptcl = 0;
	unsigned long int n_ptcl = 0;
	int n_child = 0;
	double box_size = 0;
} Box;

typedef struct Colleague{
	unsigned long int clg[27] = {0};
} Colleague;


//typedef struct Box{
//	double center[3] = {0,0,0};
//	Box * parent = NULL;
//	Box * child[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//	unsigned long int first_charge = 0;
//	int number_of_child = 0;
//	double box_size = 0;
//} Box;


int create_tree(double * x, double * y, double * z, const unsigned long int n, const unsigned int s, vector<Box> &tree, unsigned long int * list);
int create_colleague(vector<Box> &box, vector<Colleague> &clg);
int separate(Box &box1, Box &box2);

std::ostream& operator<<(std::ostream& os, Box& box);
std::ostream& operator<<(std::ostream& os, Colleague& clg);
#endif
