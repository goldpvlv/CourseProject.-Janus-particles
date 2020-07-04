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



template<typename T> struct Box : std::unique_ptr<T> {
	using std::unique_ptr<T>::unique_ptr;
	operator T &() const { return **this; }
};



class System
{
public:

	string source;
	vector <double> u;
	vector <double> grad;

	void ReadParameters();

	int layers_x;
	int layers_y;
	double Rcurv, chi_seg;
	string geometry_name;
	System(string _source);
	Geometry geo;
	void SetGeometry();
	vector<Molecule> mol;
	void ReadMolecules();
	vector<Box<BaseOptimTools>> methods;
	void ReadMethods();
	void Function();
	void Cycling();



};