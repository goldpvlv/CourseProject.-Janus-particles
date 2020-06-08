#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>

using namespace std;

const double pi = 3.14159265359;
using vector3d = vector<vector<vector<double>>>;

class Geometry
{
	public:
	void AllocateMemory(int Mx, int My);
	void Print(int Mx, int My);

	vector < vector < double > > lambda_bb;
	vector < vector < double > > lambda_bn;
	vector < vector < double > > lambda_bf;
	vector < vector < double > > lambda_fb;
	vector < vector < double > > lambda_fn;
	vector < vector < double > > lambda_ff;
	vector < vector < double > > lambda_nb;
	vector < vector < double > > lambda_nn;
	vector < vector < double > > lambda_nf;

};


class Polar : public  Geometry
{
	public:
		
	void GetValue(int Mx, int My);
	void MemoryVectors();
	void UpdateVolume();
	void UpdateSquareUp();
	void UpdateSquareSide();
	void UpdateSquareFront();
	void TranspositionPolar();

	int Mx, My;
	vector < vector < double > > volume;
	vector < double > square_up;
	vector < double > square_side;
	vector < double > square_front;
};



class Molecule
{

};


class System
{
public:
	int Mx, My, N;
	double theta, xmin, xmax, ymin, ymax;
	int M = (Mx + 2)*(My + 2);
	double teta = N * theta;

	vector < Molecule > Molecules;
	vector <double> grad;
	vector <double> direction;
	vector <double> alpha;
	vector <double> beta;
	vector <double> u;
	vector < vector < double > > A;

	void GetValue(int Mx, int My,int N,double theta, double xmin,
		double xmax, double ymin, double ymax);
	void SingularMatrix();
};



class SCF : public System {
	public:
	void MemoryVectors();
	
	double q;
	vector3d Gforw;
	vector3d Gback;
	vector < vector < double > > fi_p;
	vector < vector < double > > fi_w;
	vector < vector < double > > G;
	vector3d FindGback(Geometry& lambda_bb, Geometry& lambda_bn, Geometry& lambda_nb, Geometry& lambda_nn,
		Geometry& lambda_nf, Geometry& lambda_bf, Geometry &lambda_fb, Geometry & lambda_fn, Geometry & lambda_ff);
	vector3d FindGforw(Geometry& lambda_bb, Geometry& lambda_bn, Geometry& lambda_nb, Geometry& lambda_nn,
		Geometry& lambda_nf, Geometry& lambda_bf, Geometry &lambda_fb, Geometry & lambda_fn, Geometry & lambda_ff);
	vector <vector<double>> FindG();
	double FindQ(Polar &volume);
	vector <vector<double>> FindFiP();
	vector <vector<double>> FindFiW();
	vector <double> FindGrad();
	void MemoryVectors();

};


