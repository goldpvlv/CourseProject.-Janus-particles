#include <iostream>
#include <cmath>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>


#include "geometry.h"


using namespace std;





void Geometry::Print() {

	for (int i = 0; i < Mx; ++i) {
		for (int j = 0; j < My; ++j) {
			cout << lambda_bb[i][j] << " ";
		}
		cout << endl;
	}
}

void Geometry::AllocateMemory() {


	lambda_bb.assign(Mx + 2, vector<double>(My + 2, 0));
	lambda_bn.assign(Mx + 2, vector<double>(My + 2, 0));
	lambda_bf.assign(Mx + 2, vector<double>(My + 2, 0));
	lambda_fb.assign(Mx + 2, vector<double>(My + 2, 0));
	lambda_fn.assign(Mx + 2, vector<double>(My + 2, 0));
	lambda_ff.assign(Mx + 2, vector<double>(My + 2, 0));
	lambda_nb.assign(Mx + 2, vector<double>(My + 2, 0));
	lambda_nn.assign(Mx + 2, vector<double>(My + 2, 0));
	lambda_nf.assign(Mx + 2, vector<double>(My + 2, 0));
	volume.assign(Mx + 2, vector<double>(My + 2, 0));;
	square_side.assign(Mx + 2, vector<double>(My + 2, 0));
	square_up.assign(Mx + 2, vector<double>(My + 2, 0));
	square_front.assign(Mx + 2, vector<double>(My + 2, 0));
	square_right.assign(Mx + 2, vector<double>(My + 2, 0));;
	square_left.assign(Mx + 2, vector<double>(My + 2, 0));;

};


void Geometry::GetValue(int Mx_, int My_, double R0_) {
	Mx = Mx_;
	My = My_;
	R0 = R0_;
	double alpha = 360 / My;
	double delta_alpha = pi * alpha / 180;
}


void Geometry::PrintVolume() {
	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j) {
			cout << lambda_nn[i][j] << " ";
		}
		cout << endl;
	}
};


void Polar::Transposition() {

	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j) {
			lambda_bb[i][j] = 0.0;
			lambda_bn[i][j] = square_front[i - 1][j] / volume[i][j] * 1.0 / 6.0;
			lambda_bf[i][j] = 0.0;

			lambda_fb[i][j] = 0.0;
			lambda_fn[i][j] = square_front[i][j] / volume[i][j] * 1.0 / 6.0;
			lambda_ff[i][j] = 0.0;


			double prob = 1.0 - lambda_bn[i][j] - lambda_fn[i][j] - 1.0 / 6.0 * square_front[i][j] / volume[i][j] - 1.0 / 6.0 * square_front[i - 1][j] / volume[i][j];
			double sum = 2.0 * (square_side[i][j] + volume[i][j]);

			lambda_nb[i][j] = prob * square_side[i][j] / sum;
			lambda_nn[i][j] = prob * 2.0 * volume[i][j] / sum;
			lambda_nf[i][j] = prob * square_side[i][j] / sum;

		}
	}
}

void Polar::UpdateSquareFront() {
	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j) {
			square_front[i][j] = 0.5*My*(i *i - (i - 1) * (i - 1));
		}
	}
}

void Polar::UpdateSquareSide() {

	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j) {
			square_side[i][j] = delta_h * (i *i - (i - 1) * (i - 1));
		}
	}
}

void Polar::UpdateSquareUp() {

	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j) {
			square_up[i][j] = delta_h * i * My;
		}
	}
}

void Polar::UpdateVolume() {

	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j) {
			volume[i][j] = 0.5*My*(i *i - (i - 1) * (i - 1))*delta_h;
		}
	}
}

void Torus::UpdateSquareFront() {


	for (int i = 1; i < Mx + 1; i++)
		for (int j = 1; j < My + 1; j++) {
			square_front[i][j] = i * delta_alpha * 2.0 * pi*(R0 + i * cos((2 * j - 1)*delta_alpha) / 2.0);
		}

}

void Torus::UpdateVolume() {
	vector <double> rho(Mx + 1, 0);

	for (int i = 1; i < Mx + 1; ++i) {
		rho[i] = 0.67*(2 * i - 1);
	}

	for (int i = 1; i < Mx + 1; i++) {
		for (int j = 1; j < My + 1; j++) {
			volume[i][j] = square_up[i][j] * pi * 2.0*(R0 + rho[i] * cos(((2 * j - 1) / 2.0)*delta_alpha));
		}
	}

};

void Torus::UpdateSquareUp() {

	for (int i = 1; i < Mx + 1; i++)
		for (int j = 1; j < My + 1; j++) {
			square_up[i][j] = delta_alpha / 2.0*(i *i - (i - 1) * (i - 1));
		}
};

void Torus::UpdateSquareRight() {

	for (int i = 1; i < Mx + 1; i++)
		for (int j = 1; j < My + 1; j++) {
			square_right[i][j] = (i - (i - 1)) * 2.0 * pi*(R0 + ((2 * i - 1) / 2.0 * cos(j*delta_alpha)));
		}

};

void Torus::UpdateSquareLeft() {

	for (int i = 1; i < Mx + 1; i++)
		for (int j = 1; j < My + 1; j++) {
			square_left[i][j] = (i - (i - 1)) * 2.0 * pi*(R0 + ((2 * i - 1) / 2.0*cos((j - 1)*delta_alpha)));
		}

};

void Torus::Transposition() {

	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j) {
			lambda_bb[i][j] = 0;
			lambda_bn[i][j] = square_front[i - 1][j] / volume[i][j] * 1.0 / 6.0;
			lambda_bf[i][j] = 0.0;

			lambda_fb[i][j] = 0.0;
			lambda_fn[i][j] = square_front[i][j] / volume[i][j] * 1.0 / 6.0;
			lambda_ff[i][j] = 0.0;


			double prob = 1.0 - 1.0 / 6.0 * square_front[i - 1][j] / volume[i][j] - 1.0 / 6.0 * square_front[i - 1][j] / volume[i][j];
			double sum = square_left[i][j] + square_right[i][j] + 2 * volume[i][j];

			lambda_nb[i][j] = prob * square_left[i][j] / sum;
			lambda_nn[i][j] = prob * 2.0 * volume[i][j] / sum;
			lambda_nf[i][j] = prob * square_right[i][j] / sum;

		}
	}

};