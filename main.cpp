#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include "classes.h"

using namespace std;

int main()
{
	double  theta=0.01, del, nu, max_step, xmin, xmax, ymin, ymax, q = 0, R_g = 0, k = 0;
	int N=10, Mx=10, My=10;
	string key="polar";

	//freopen("Text.txt", "r", stdin);

	//cin >> N >> Mx >> My >> theta >> nu >> del >> max_step >> key >> xmin >> xmax >> ymin >> ymax;
	 
	int M = (Mx + 2)*(My + 2);
	double teta = N * theta;

	if (key == "polar") {

		Polar polar;

		polar.GetValue(Mx, My);
		polar.AllocateMemory(Mx, My);
		polar.MemoryVectors();
		polar.UpdateVolume();
		polar.UpdateSquareUp();
		polar.UpdateSquareSide();
		polar.UpdateSquareFront();
		polar.TranspositionPolar();
		polar.Print(Mx,My);

	}
	else if (key == "torus") {

	}
	else if (key == "sphere") {

	}

	SCF grad;

	grad.GetValue(Mx, My, N, theta, xmin, xmax, ymin, ymax);
	grad.SingularMatrix();
	grad.MemoryVectors();
	grad.FindG();
	grad.FindGforw(Geometry& lambda_bb, Geometry& lambda_bn, Geometry& lambda_nb, Geometry& lambda_nn,
		Geometry& lambda_nf, Geometry& lambda_bf, Geometry &lambda_fb, Geometry & lambda_fn, Geometry & lambda_ff);
	grad.FindGback(Geometry& lambda_bb, Geometry& lambda_bn, Geometry& lambda_nb, Geometry& lambda_nn,
		Geometry& lambda_nf, Geometry& lambda_bf, Geometry &lambda_fb, Geometry & lambda_fn, Geometry & lambda_ff);
	grad.FindQ(Polar &volume);
	grad.FindFiP();
	grad.FindFiW();
	grad.FindGrad();



	

	return 0;
}