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

	virtual double SetGradFirst(vector <double> grad)  { double q = 0.0; return q; };
	virtual double SetGradRegular(vector <double> grad) { double q = 0.0; return q; };
	virtual void UpdateX(vector <double> &x, vector <double> grad) {};

	void SetParameters(int Mx, int My, vector <double> grad, vector <double> u);

	virtual ~BaseOptimTools() {};
};


class Gradient : public BaseOptimTools
{
public:
	double SetGradFirst(vector <double> grad) override;
	double SetGradRegular(vector <double> grad) override;
	void UpdateX(vector <double> &x, vector <double> grad) override;

	Gradient(int _num_iter, double _tolerance, double _nu);

};


class DFP : public BaseOptimTools
{
public:

	double SetGradFirst(vector <double> grad) override;
	double SetGradRegular(vector <double> grad) override;
	void UpdateX(vector <double> &x, vector <double> grad) override;

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