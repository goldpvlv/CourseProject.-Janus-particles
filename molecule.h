#pragma once


#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
using vector3d = vector<vector<vector<double>>>;


class Molecule {
public:

	double theta, sigma, chi;
	int xmin, xmax, ymin, ymax, num_generation, ns, q;
	
	void MyNewMethod();
	void AllocateMemory();

	double q;
	vector3d Gback;
	vector3d Gforw;
	vector <double> u;
	vector <vector <double>> G;
	vector <double> FiSide;
	vector<vector<double>> fi_p;
	vector<vector<double>> fi_w;
	vector<vector<double>> fi;
	vector<vector<double>> multipliers;



};