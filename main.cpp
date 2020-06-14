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
	double  theta, del, nu, max_step, xmin, xmax, ymin, ymax, q =0, R_g = 0, k = 0;
	int N, Mx, My;
	string key="polar";

	freopen("Text.txt", "r", stdin);

	cin >> N >> Mx >> My >> theta >> nu >> del >> max_step >> key >> xmin >> xmax >> ymin >> ymax;
	int M = (Mx + 2)*(My + 2);
	double teta = N * theta;
	Polar polar;
	if (key == "polar") {

		

		polar.GetValue(Mx, My);
		polar.AllocateMemory(Mx, My);
		polar.MemoryVectors();
		polar.UpdateVolume();
		polar.UpdateSquareUp();
		polar.UpdateSquareSide();
		polar.UpdateSquareFront();
		polar.PrintVolume();
		cout << endl;
		polar.TranspositionPolar(polar.lambda_bb, polar.lambda_bn, polar.lambda_nb, polar.lambda_nn,
			polar.lambda_nf, polar.lambda_bf, polar.lambda_fb, polar.lambda_fn, polar.lambda_ff);
		//polar.Print(Mx,My);
		cout << endl;

	}
	else if (key == "torus") {

	}
	else if (key == "sphere") {

	}

	SCF gradient;

	gradient.GetValue(Mx, My, N, M, theta, xmin, xmax, ymin, ymax);
	gradient.MemoryVectors();
	gradient.SingularMatrix();
	gradient.FindG();
	gradient.FindGforw(polar.lambda_bb, polar.lambda_bn, polar.lambda_nb, polar.lambda_nn,
		polar.lambda_nf, polar.lambda_bf, polar.lambda_fb, polar.lambda_fn, polar.lambda_ff);
	gradient.FindGback(polar.lambda_bb, polar.lambda_bn, polar.lambda_nb, polar.lambda_nn,
		polar.lambda_nf, polar.lambda_bf, polar.lambda_fb, polar.lambda_fn, polar.lambda_ff);
	gradient.FindQ(polar.volume);
	gradient.FindFiP();
	gradient.FindFiW();
	gradient.FindGrad();
	//gradient.Print();



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


		gradient.GetValue(Mx, My, N, M, theta, xmin, xmax, ymin, ymax);
		gradient.MemoryVectors();
		gradient.SingularMatrix();
		gradient.FindG();
		gradient.FindGforw(polar.lambda_bb, polar.lambda_bn, polar.lambda_nb, polar.lambda_nn,
			polar.lambda_nf, polar.lambda_bf, polar.lambda_fb, polar.lambda_fn, polar.lambda_ff);
		gradient.FindGback(polar.lambda_bb, polar.lambda_bn, polar.lambda_nb, polar.lambda_nn,
			polar.lambda_nf, polar.lambda_bf, polar.lambda_fb, polar.lambda_fn, polar.lambda_ff);
		gradient.FindQ(polar.volume);
		gradient.FindFiP();
		gradient.FindFiW();
		gradient.FindGrad();

		for (int i = 0; i < M; i++)
			gradient.beta[i] = gradient.grad[i] - gradient.beta[i];


		gradient.Formula();

		k++;

	} while ((gradient.length_of_grad > del) && (k < max_step));


	FILE * txt = fopen("fi_p_out.txt", "w");
	fprintf(txt, "%5d ", (int)My);
	for (int j = 1; j < My + 1; j++) {
		fprintf(txt, " %11d", j);
	}
	fprintf(txt, "\n");
	for (int i = 1; i < Mx + 1; i++) {
		fprintf(txt, "%5d ", i);
		for (int j = 1; j < My + 1; j++) {
			fprintf(txt, " %8.5e", gradient.fi_p[i][j]);
		}
		fprintf(txt, "\n");
	}
	fclose(txt);


	_getch();

	return 0;
}