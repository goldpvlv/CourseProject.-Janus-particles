#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include <conio.h>

#include "molecule.h"


using namespace std;

void Molecule::MyNewMethod() {

	int num_atoms = 1 + ns * (pow(2, num_generation + 1) - 1);
	double amount_of_substance = sigma * num_atoms;

};


