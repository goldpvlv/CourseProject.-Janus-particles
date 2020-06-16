#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include "classes.h"
#include <conio.h>

using namespace std;

int main()
{
	double  theta, tolerance, nu, num_iter, xmin, xmax, ymin, ymax, q =0, R_g = 0, step = 0;
	int N, layers_x, layers_y;
	string key="polar";

	freopen("input.txt", "r", stdin);

	cin >> N >> layers_x >> layers_y >> theta >> nu >> tolerance >> num_iter >> key >> xmin >> xmax >> ymin >> ymax;
	int M = (layers_x + 2)*(layers_y + 2);
	double teta = N * theta;

	Geometry geo;

	if (key == "polar") {
		Polar polar;
		polar.GetValue(layers_x, layers_y, N);
		polar.AllocateMemory();
		polar.UpdateVolume();
		polar.UpdateSquareUp();
		polar.UpdateSquareSide();
		polar.UpdateSquareFront();
		polar.PrintVolume();
		cout << endl;
		polar.Transposition(polar.lambda_bb, polar.lambda_bn, polar.lambda_nb, polar.lambda_nn,
			polar.lambda_nf, polar.lambda_bf, polar.lambda_fb, polar.lambda_fn, polar.lambda_ff, polar.volume, 
			polar.square_up, polar.square_side, polar.square_front);
		polar.Print();


		cout << endl;

	}
	else if (key == "torus") {

		Torus torus;



	}
	else if (key == "sphere") {

	}

	Molecule molecule;
	System gradient;

	molecule.GetValue(layers_x, layers_y, N, theta, xmin, xmax, ymin, ymax);
	molecule.MemoryVectors();
	gradient.SingularMatrix();
	gradient.FindG();	
	gradient.FindGforw(geo.lambda_bb, geo.lambda_bn,geo.lambda_nb, geo.lambda_nn,
		geo.lambda_nf, geo.lambda_bf, geo.lambda_fb, geo.lambda_fn, geo.lambda_ff);
	gradient.FindGback(geo.lambda_bb, geo.lambda_bn, geo.lambda_nb, geo.lambda_nn,
		geo.lambda_nf, geo.lambda_bf, geo.lambda_fb, geo.lambda_fn, geo.lambda_ff);
	gradient.FindQ(geo.volume);
	gradient.FindFiP();
	gradient.FindFiW();
	gradient.SingularMatrix();
	gradient.FindGrad(molecule.fi_p, molecule.fi_w);
	gradient.Print();



	do {

		gradient.FindDirection();

		for (int i = 0; i < M; i++)                        
			gradient.alpha[i] = gradient.u[i];

		for (int i = 0; i < M; i++)                        
			gradient.u[i] = gradient.u[i] + nu * gradient.direction[i];


		for (int i = 0; i < M; i++)                        
			gradient.alpha[i] = gradient.u[i] - gradient.alpha[i];

		for (int i = 0; i < M; i++)                      
			gradient.beta[i] = gradient.grad[i];


		gradient.GetValue(layers_x, layers_y, N, M, theta, xmin, xmax, ymin, ymax);
		gradient.MemoryVectors();
		gradient.SingularMatrix();
		gradient.FindG();
		gradient.FindGforw(geo.lambda_bb, geo.lambda_bn, geo.lambda_nb, geo.lambda_nn,
			geo.lambda_nf, geo.lambda_bf, geo.lambda_fb, geo.lambda_fn, geo.lambda_ff);
		gradient.FindGback(geo.lambda_bb, geo.lambda_bn, geo.lambda_nb, geo.lambda_nn,
			geo.lambda_nf, geo.lambda_bf, geo.lambda_fb, geo.lambda_fn, geo.lambda_ff);
		gradient.FindQ(geo.volume);
		gradient.FindFiP();
		gradient.FindFiW();
		gradient.FindGrad();

		for (int i = 0; i < M; i++)
			gradient.beta[i] = gradient.grad[i] - gradient.beta[i];


		gradient.Formula();

		step++;

	} while ((gradient.length_of_grad > tolerance) && (step < num_iter));


	FILE * txt = fopen("fi_p_out.txt", "w");
	fprintf(txt, "%5d ", (int)layers_x);
	for (int j = 1; j < layers_x + 1; j++) {
		fprintf(txt, " %11d", j);
	}
	fprintf(txt, "\n");
	for (int i = 1; i < layers_x + 1; i++) {
		fprintf(txt, "%5d ", i);
		for (int j = 1; j < layers_y + 1; j++) {
			fprintf(txt, " %8.5e", gradient.fi_p[i][j]);
		}
		fprintf(txt, "\n");
	}
	fclose(txt);


	_getch();

	return 0;
}