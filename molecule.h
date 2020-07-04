#pragma once


#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


class Molecule {
public:

	double theta, sigma, chi;
	int xmin, xmax, ymin, ymax, num_generation, ns, q;
	
	void MyNewMethod();

};