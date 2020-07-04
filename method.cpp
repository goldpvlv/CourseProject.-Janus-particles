#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include <conio.h>

#include "method.h"


using namespace std;

void BaseOptimTools::SetParameters(double t, int ns, double eta) {

	M = 3 * (Mx + 2)*(My + 2);
};

Gradient::Gradient(int _num_iter, double _tolerance, double _nu) {
	name = "Gradient";
	num_iter = _num_iter;
	tolerance = _tolerance;
	nu = _nu;
}

double Gradient::SetGradFirst(vector <double> grad) {

	double lenght=0.0;

	for (int i = 0; i < M; i++)
		lenght += grad[i]*grad[i];
	return lenght;
};

double Gradient::SetGradRegular(vector <double> grad) {

	double lenght = 0.0;

	for (int i = 0; i < M; i++)
		lenght += grad[i] * grad[i];
	return lenght;
};

void Gradient::UpdateX(vector <double> u, vector <double> grad) {
	for (int i = 0; i < grad.size(); i++)
		u[i] = u[i] - 0.001 * grad[i];
};


DFP::DFP(int _num_iter, double _tolerance, double _nu) {
	name = "DFP";
	num_iter = _num_iter;
	tolerance = _tolerance;
	nu = _nu;
}

double DFP::SetGradFirst(vector <double> grad) {

	double lenght = 0.0;

	for (int i = 0; i < M; i++)
		lenght += grad[i] * grad[i];
	return lenght;

};


double DFP::SetGradRegular(vector <double> grad) {

	for (int i = 0; i < grad.size(); i++)
		beta[i] = grad[i] - beta[i];
	A = Formula();

	double lenght = 0.0;

	for (int i = 0; i < M; i++)
		lenght += grad[i] * grad[i];
	return lenght;

};



void DFP::UpdateX(vector <double> u, vector <double> grad) {

	direction = FindDirection(grad);
	for (int i = 0; i < grad.size(); i++)
		alpha[i] = u[i];
	for (int i = 0; i < grad.size(); i++)
		u[i] = u[i] + nu * direction[i];
	for (int i = 0; i < grad.size(); i++)
		alpha[i] = u[i] - alpha[i];

	for (int i = 0; i < grad.size(); i++)
		beta[i] = grad[i];


}

void DFP::SingularMatrix() {
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			if (i == j)
				A[i][j] = 1;
			else A[i][j] = 0;
		}
	}
};

vector<double> DFP::FindDirection(vector<double>grad) {

	vector <double> direction(M, 0);

	direction = Multiply_matrix_by_vector(A, grad);

	for (int i = 0; i < M; i++)
		direction[i] = -1 * direction[i];
	return direction;

};


vector < vector < double > > DFP::Formula() {
	double number1 = 0, number2 = 0;
	vector < vector < double > > matrix1(M, vector<double>(M, 0));
	vector < vector < double > > matrix2(M, vector<double>(M, 0));
	vector < vector < double > > B(M, vector<double>(M, 0));
	vector < vector < double > > C(M, vector<double>(M, 0));

	vector <double> a(M, 0);
	vector <double> b(M, 0);


	matrix1 = Multiply_vectors(alpha, alpha);
	number1 = Find_number(alpha, beta);
	B = Division_matrix_on_number(matrix1, number1);
	a = Multiply_matrix_by_vector(A, beta);
	matrix2 = Multiply_vectors(a, beta);
	matrix2 = Multiply_matrixes(A, matrix1);
	b = Multiply_matrix_by_vector(A, beta);
	number2 = Find_number(b, beta);
	C = Division_matrix_on_number(matrix2, number2);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++)
			A[i][j] = A[i][j] + B[i][j] - C[i][j];
	}

	return A;
};

vector < vector < double > > DFP::Multiply_matrixes(vector < vector < double > > A, vector < vector < double > > U) {

	vector < vector < double > > c(M, vector<double>(M, 0));

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			c[i][j] = 0;
			for (int t = 0; t < M; t++)
				c[i][j] += A[i][t] * U[t][j];
		}
	}
	return c;
};


vector < vector < double > > DFP::Multiply_vectors(vector<double> a, vector<double> b) {

	vector < vector < double > > matrix(M, vector<double>(M, 0));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++)
			matrix[i][j] = a[i] * b[j];
	}
	return matrix;
};

double DFP::Find_number(vector<double> a, vector<double> b) {
	double number = 0;
	for (int i = 0; i < M; i++)
		number += a[i] * b[i];
	return number;
};


vector < double > DFP::Multiply_matrix_by_vector(vector<vector<double>>A, vector<double>grad) {
	vector <double> a(M, 0);
	for (int i = 0; i < M; i++) {
		a[i] = 0;
		for (int j = 0; j < M; j++) {
			a[i] += A[i][j] * grad[j];
		}
	}
	return a;
};

vector < vector < double > > DFP::Division_matrix_on_number(vector<vector<double>>A, double b) {
	vector < vector < double > > matrix1(M, vector<double>(M, 0));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++)
			matrix1[i][j] = A[i][j] / b;
	}
	return matrix1;
};
