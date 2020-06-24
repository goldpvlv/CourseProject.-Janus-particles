#pragma once


#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

const double pi = 3.14159265359;
using vector3d = vector<vector<vector<double>>>;

template<typename T> struct Box : std::unique_ptr<T> {
	using std::unique_ptr<T>::unique_ptr;
	operator T &() const { return **this; }
};


class BaseOptimTools
{
public:
	int Mx, My, R0, xmax, xmin, ymax, ymin, theta;
	int num_iter;
	double tolerance, q, length_of_grad, nu;
	double step;
	double teta = R0 * theta;
	int M = 3*(Mx + 2)*(My + 2);

	string name = "OptimTools";
	BaseOptimTools(int num_iter = 0, double tolerance = 0, double step = 0) noexcept : num_iter(num_iter), tolerance(tolerance), step(step) { };
	BaseOptimTools(int);
	void GetValue(int Mx_, int My_, int R0);
	void AllocateMemory();
	void Print();

	virtual void SingularMatrix() const;
	virtual vector <vector <double>> FindG() const;
	virtual vector3d FindGforw() const;
	virtual vector3d FindGback() const;
	virtual vector <vector<double>> FindFiW() const;
	virtual vector <vector<double>> FindFiP() const;
	virtual vector <double> FindGrad() const;
	virtual double FindQ()  const;
	virtual double LengthOfGrad()const;
	virtual vector < vector < double > > Formula() const;
	vector < double >  Multiply_matrix_by_vector(vector < vector < double > >, vector<double>);

	vector<double> FindDirection(vector<vector<double>>, vector<double>);
	vector <double> Function();

	void Сycling();


	//friend BaseOptimTools operator/(int a, const BaseOptimTools &matrix); //division_matrix_by_number
	friend BaseOptimTools operator*(const BaseOptimTools& matrix);// multiply 2 matrixes


	vector3d Gforw;
	vector3d Gback;
	vector < vector < double > >fi_p;
	vector < vector < double > >fi_w;
	vector < vector < double > > G;
	vector < vector < double > > A;
	vector < vector < double > > arr;
	vector<double> p;
	vector <double> grad;
	vector <double> direction;
	vector <double> alpha;
	vector <double> beta;
	vector <double> u;
};


class Gradient : public BaseOptimTools
{
public:
	Gradient(double t); // конструктор
	vector <vector <double>> FindG();
	vector3d FindGforw();
	vector3d FindGback();
	vector <vector<double>> FindFiW();
	vector <vector<double>> FindFiP();
	vector <double> FindGrad();
	double FindQ();
	double LengthOfGrad();

};


class DFP : public BaseOptimTools
{
public:


	DFP(double t);
	void SingularMatrix();
	vector < vector < double > > Formula();
	vector < vector < double > > Multiply_matrixes(vector < vector < double > >, vector < vector < double > >);
	vector < vector < double > > Multiply_vectors(vector<double>, vector<double>);
	double Find_number(vector<double>, vector<double>);
	vector < double >  Multiply_matrix_by_vector(vector < vector < double > >, vector<double>);
	vector < vector < double > > Division_matrix_on_number(vector < vector < double >>, double);


};



class Geometry {
public:
	mutable int Mx, My;
	double R0;
	void GetValue(int Mx_, int My_, double R0_);
	double alpha = 360 / My;
	double delta_alpha = pi * alpha / 180;
	//  

	vector < vector < double > >lambda_bb;
	vector < vector < double > >lambda_bn;
	vector < vector < double > >lambda_bf;
	vector < vector < double > >lambda_fb;
	vector < vector < double > >lambda_fn;
	vector < vector < double > >lambda_ff;
	vector < vector < double > >lambda_nb;
	vector < vector < double > >lambda_nn;
	vector < vector < double > >lambda_nf;
	vector < vector < double > >volume;
	vector < vector < double > >square_side;
	vector < vector < double > >square_up;
	vector < vector < double > >square_front;
	vector < vector < double > >square_right;
	vector < vector < double > >square_left;
	virtual void UpdateVolume() const; // только объявляем: реализация в дочерних классах
	virtual void UpdateSquareFront() const;
	virtual void UpdateSquareUp() const;
	virtual void UpdateSquareLeft() const;
	virtual void UpdateSquareRight() const;
	virtual void UpdateSquareSide() const;
	virtual void Transposition() const;
	void PrintVolume();
	void Print();
	//
	void AllocateMemory();
};



class Polar : public Geometry {
public:
	void UpdateVolume();
	void UpdateSquareFront();
	void UpdateSquareUp();
	void UpdateSquareSide();
	void Transposition();

};

class Torus : public Geometry {
public:
	void UpdateVolume();
	void UpdateSquareFront();
	void UpdateSquareUp();
	void UpdateSquareLeft();
	void UpdateSquareRight();
	void Transposition();
};



class Molecule {
public:


	int num_generation, line_item, density_of_molecules, statical_sum;
	int num_atoms = 1 + line_item * (pow(2, num_generation + 1) - 1);
	double amount_of_substance = density_of_molecules * num_atoms;



};


class System
{
public:
	// название входного файла:
	string source;
	//
	void ReadParameters();

	int layers_x;
	int layers_y;
	double Rcurv;
	//// про Геометрию
	string geometry_name;
	Geometry geo;
	void SetGeometry();
	// молекулы 
	vector<Molecule> mol;
	void ReadMolecules();
	// методы сходимости 
	vector<Box<BaseOptimTools>> methods;
	void ReadMethods();
	// конструктор
	System(string _source);

	void GetMethods();

};


