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
	void GetValue(int Mx, int My, int N);
	void AllocateMemory();
	void Print();


	vector < vector < double > > lambda_bb;
	vector < vector < double > > lambda_bn;
	vector < vector < double > > lambda_bf;
	vector < vector < double > > lambda_fb;
	vector < vector < double > > lambda_fn;
	vector < vector < double > > lambda_ff;
	vector < vector < double > > lambda_nb;
	vector < vector < double > > lambda_nn;
	vector < vector < double > > lambda_nf;
	int Mx, My, N;
	vector < vector < double > > volume;
	vector < vector < double > > square_side;
	vector < vector < double > > square_up;
	vector < vector < double > > square_right;
	vector < vector < double > > square_left;
	vector < vector < double > > square_front;


};


class Polar : public  Geometry
{
	public:
	
	void UpdateVolume();
	void UpdateSquareUp();
	void UpdateSquareSide();
	void UpdateSquareFront();
	void Transposition(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn,
		vector<vector<double>> & lambda_nb, vector<vector<double>> & lambda_nn, vector<vector<double>> & lambda_nf,
		vector<vector<double>> & lambda_bf, vector<vector<double>> &lambda_fb, vector<vector<double>> & lambda_fn, 
		vector<vector<double>> & lambda_ff, vector<vector<double>> & volume,vector<vector<double>> & square_up, 
		vector<vector<double>> & square_side, vector<vector<double>> & square_front);

	void PrintVolume();

	
};


class Torus :public Geometry {
public:

	double alpha = 360 / My;
	double sum1 = 0, sum2 = 0;
	double delta_alpha = pi * alpha / 180;

	void UpdateVolume();
	void UpdateSquareUp();
	void UpdateSquareRight();
	void UpdateSquareLeft();
	void UpdateSquareFront();
	void Transposition(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn,
		vector<vector<double>> & lambda_nb, vector<vector<double>> & lambda_nn, vector<vector<double>> & lambda_nf,
		vector<vector<double>> & lambda_bf, vector<vector<double>> &lambda_fb, vector<vector<double>> & lambda_fn,
		vector<vector<double>> & lambda_ff, vector<vector<double>> & volume, vector<vector<double>> & square_up,
		vector<vector<double>> & square_left, vector<vector<double>> & square_right,  vector<vector<double>> & square_front);


};

class Sperical :public Geometry {

};

class Molecule {


	int num_generation, line_item, density_of_molecules, statical_sum ;
	int num_atoms = 1 + line_item*(pow(2, num_generation + 1) - 1);
	double amount_of_substance = density_of_molecules * num_atoms;



	void GetValue(int Mx, int My, int N, double theta, double xmin,	double xmax, double ymin, double ymax);
	void MemoryVectors();
	vector <vector<double>> FindG(vector <double> & u);
	vector <vector<double>> FindFiP();
	vector <vector<double>> FindFiW();
	vector3d FindGback(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn,
		vector<vector<double>> & lambda_nb, vector<vector<double>> & lambda_nn, vector<vector<double>> & lambda_nf,
		vector<vector<double>> & lambda_bf, vector<vector<double>> &lambda_fb, vector<vector<double>> & lambda_fn, vector<vector<double>> & lambda_ff);
	vector3d FindGforw(vector<vector<double>> & lambda_bb, vector<vector<double>> &lambda_bn, vector<vector<double>> &lambda_nb,
		vector<vector<double>> &lambda_nn, vector<vector<double>> &lambda_nf, vector<vector<double>> &lambda_bf,
		vector<vector<double>> &lambda_fb, vector<vector<double>> & lambda_fn, vector<vector<double>> & lambda_ff);
	double FindQ(vector<vector<double>> & volume);




	vector3d Gforw;
	vector3d Gback;
	vector < vector < double > > fi_p;
	vector < vector < double > > fi_w;
	vector < vector < double > > G;
	
	int Mx, My, N, xmax, xmin, ymax, ymin, theta;
	double q;
	int M = (Mx + 2)*(My + 2);

};


class System
{
public:

	void GetValue(int Mx, int My);
	void MemoryVectors();
	void SingularMatrix();
	vector <double> FindGrad(vector < vector < double > > &fi_p, vector < vector < double > > &fi_w);
	vector < double >  FindDirection();
	void LengthOfGrad();
	vector < vector < double > > Formula();
	void Print();

	vector < Molecule > Molecules;
	vector <double> u;
	vector <double> grad;
	vector <double> direction;
	vector <double> alpha;
	vector <double> beta;
	vector < vector < double > > A;
	int Mx, My, M=(Mx*2)+(My*2);
	double length_of_grad;
};




