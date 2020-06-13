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
	double  theta=0.01, del, nu, max_step, xmin=1, xmax=1, ymin=1, ymax=1, q = 0, R_g = 0, k = 0;
	int N=10, Mx=10, My=10, M=(Mx+2)*(My+2);
	string key="polar";

	//freopen("Text.txt", "r", stdin);

	//cin >> N >> Mx >> My >> theta >> nu >> del >> max_step >> key >> xmin >> xmax >> ymin >> ymax;
	 
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
		polar.Print(Mx,My);
		cout << endl;

	}
	else if (key == "torus") {

	}
	else if (key == "sphere") {

	}

	SCF grad;

	grad.GetValue(Mx, My, N, M, theta, xmin, xmax, ymin, ymax);
	grad.MemoryVectors();
	grad.SingularMatrix();
	grad.FindG();
	grad.FindGforw(polar.lambda_bb, polar.lambda_bn, polar.lambda_nb, polar.lambda_nn,
		polar.lambda_nf, polar.lambda_bf, polar.lambda_fb, polar.lambda_fn, polar.lambda_ff);
	grad.FindGback(polar.lambda_bb, polar.lambda_bn, polar.lambda_nb, polar.lambda_nn,
		polar.lambda_nf, polar.lambda_bf, polar.lambda_fb, polar.lambda_fn, polar.lambda_ff);
	grad.FindQ(polar.volume);
	grad.FindFiP();
	grad.FindFiW();
	grad.FindGrad();
	grad.Print();



	

	return 0;
}