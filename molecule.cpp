#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include <conio.h>

#include "molecule.h"



using namespace std;

void Molecule::SetParameters() {

	num_atoms = 1 + ns * (pow(2, num_generation + 1) - 1);
	theta = sigma * num_atoms;

};


void Molecule::AllocateMemory(int layers_x_, int layers_y_, int M) {

	layers_x = layers_x_;
	layers_y = layers_y_;

	Gback.assign(layers_x + 2, vector <vector<double>>(layers_y + 2, vector<double>(num_atoms + 1, 0)));
	Gforw.assign(layers_x + 2, vector <vector<double>>(layers_y + 2, vector<double>(num_atoms + 1, 0)));;
	u.resize(M+2,0);
	G.assign(layers_x + 2, vector<double>(layers_y + 2, 0));
	fi_side.resize(M+2,0);
	fi_p.assign(layers_x + 2, vector<double>(layers_y + 2, 0));
	multipliers.assign(layers_x + 1, vector <double>(layers_y + 1, 0));

};

void Molecule::FindG() {


	for (int i = 0; i < layers_x + 1; ++i) {
		for (int j = 0; j < layers_y + 1; ++j) {
			G[i][j] = exp(-u[i*(layers_y + 2) + j]);
		}
	}

};

void Molecule::FindGforw(Geometry geo) {

	for (int j = 1; j < num_atoms; ++j) {
		for (int k = 1; k < layers_x + 1; ++k) {
			for (int i = 1; i < layers_y + 1; ++i) {
				
				Gforw[i][k][j] = G[i][k] * (geo.lambda_bb[i][k] * Gforw[i - 1][k - 1][j - 1] + geo.lambda_bn[i][k] * Gforw[i - 1][k][j - 1] + geo.lambda_bf[i][k] * Gforw[i - 1][k + 1][j - 1] +
					geo.lambda_nb[i][k] * Gforw[i][k - 1][j - 1] + geo.lambda_nn[i][k] * Gforw[i][k][j - 1] + geo.lambda_nf[i][k] * Gforw[i][k + 1][j - 1] +
					geo.lambda_fb[i][k] * Gforw[i + 1][k - 1][j - 1] + geo.lambda_fn[i][k] * Gforw[i + 1][k + 1][j - 1] + geo.lambda_ff[i][k] * Gforw[i + 1][k + 1][j - 1]);

				Gforw[i][layers_x + 1][j] = Gforw[i][1][j];
				Gforw[i][0][j] = Gforw[i][layers_x][j];
			}

			Gforw[layers_y + 1][k][j] = Gforw[layers_y][k][j];
			Gforw[0][k][j] = 0;
		}
	}

};


void Molecule::FindGback(Geometry geo) {


	for (int j = num_atoms - 1; j > 0; --j) {
		for (int k = 1; k < layers_x + 1; ++k) {
			for (int i = 1; i < layers_y + 1; ++i) {
				Gback[i][k][j - 1] = G[i][k] * (geo.lambda_bb[i][k] * Gback[i - 1][k - 1][j] + geo.lambda_bn[i][k] * Gback[i - 1][k][j] + geo.lambda_bf[i][k] * Gback[i - 1][k + 1][j] +
					geo.lambda_nb[i][k] * Gback[i][k - 1][j] + geo.lambda_nn[i][k] * Gback[i][k][j] + geo.lambda_nf[i][k] * Gback[i][k + 1][j] +
					geo.lambda_fb[i][k] * Gback[i + 1][k - 1][j] + geo.lambda_fn[i][k] * Gback[i + 1][k][j] + geo.lambda_ff[i][k] * Gback[i + 1][k + 1][j]);

				Gback[i][layers_x + 1][j - 1] = Gback[i][1][j - 1];
				Gback[i][0][j - 1] = Gback[i][layers_x][j - 1];

			}

			Gback[layers_y + 1][k][j - 1] = Gback[layers_y][k][j - 1];
			Gback[0][k][j - 1] = 0;

		}
	}
};



void Molecule::FindQ(Geometry geo) {

	double sum = 0;
	q = 0;
	for (int i = 1; i < layers_x + 1; ++i) {
		for (int k = 1; k < layers_y + 1; ++k) {
			sum = 0;
			for (int j = 0; j < num_atoms; ++j)
				sum += (Gforw[i][k][j] * Gback[i][k][j]) / G[i][k];
			q += sum * geo.volume[i][k];
		}
	}

};