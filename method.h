#pragma once


#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


class BaseOptimTools
{
public:
	int Mx, My, M;
	int num_iter;
	double tolerance, length_of_grad, nu;

	string name = "OptimTools";
	vector <double> direction;

	virtual double SetGradFirst(vector <double> grad) const;
	virtual double SetGradRegular(vector <double> grad) const;
	virtual vector <double> UpdateX(vector <double> x, vector <double> grad) const;

	void SetParameters(double t, int ns, double eta);

	virtual ~BaseOptimTools() {};
};


class Gradient : public BaseOptimTools
{
public:
	double SetGradFirst(vector <double> grad);
	double SetGradRegular(vector <double> grad);
	void UpdateX(vector <double> x, vector <double> grad);

	Gradient(int _num_iter, double _tolerance, double _nu);

};


class DFP : public BaseOptimTools
{
public:

	double SetGradFirst(vector <double> grad);
	double SetGradRegular(vector <double> grad);
	void UpdateX(vector <double> x, vector <double> grad);

	DFP(int _num_iter, double _tolerance, double _nu);

	vector <double> alpha;
	vector <double> beta;
	vector < vector < double > > A;

private:
	vector < vector < double > > Formula();
	void SingularMatrix();
	vector<double> FindDirection(vector<double>grad);
	vector < vector < double > > Multiply_matrixes(vector < vector < double > >, vector < vector < double > >);
	vector < vector < double > > Multiply_vectors(vector<double>, vector<double>);
	double Find_number(vector<double>, vector<double>);
	vector < double >  Multiply_matrix_by_vector(vector < vector < double > >, vector<double>);
	vector < vector < double > > Division_matrix_on_number(vector < vector < double >>, double);
};