#pragma once


#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <stdio.h>


#include "geometry.h"


using namespace std;
using vector3d = vector<vector<vector<double>>>;



class Molecule {


	friend class System;


public:

	double theta, sigma, chi;
	int xmin, xmax, ymin, ymax, num_generation, num_atoms, ns;
	
	void SetParameters();
	void AllocateMemory(int layers_x_, int layers_y_, int M);

	double q;
	vector3d Gback;
	vector3d Gforw;
	vector <double> u;
	vector <vector <double>> G;
	vector <double> fi_side;
	vector<vector<double>> fi_p;
	vector<vector<double>> multipliers;

	
	void FindG();
	void FindGforw(Geometry geo);
	void FindGback(Geometry geo);
	void FindQ(Geometry geo);

};