#include <iostream>
#include <cmath>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include "classes.h"


using namespace std;


	void Polar::TranspositionPolar() {
		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j) {
				lambda_bb[i][j] = 0.0;
				lambda_bn[i][j] = square_front[i - 1] / volume[i][j] * 1.0 / 6.0;
				lambda_bf[i][j] = 0.0;

				lambda_fb[i][j] = 0.0;
				lambda_fn[i][j] = square_front[i] / volume[i][j] * 1.0 / 6.0;
				lambda_ff[i][j] = 0.0;


				double prob = 1.0 - lambda_bn[i][j] - lambda_fn[i][j] - 1.0 / 6.0 * square_front[i] / volume[i][j] - 1.0 / 6.0 * square_front[i - 1] / volume[i][j];
				double sum = 2.0 * (square_side[i] + volume[i][j]);

				lambda_nb[i][j] = prob * square_side[i] / sum;
				lambda_nn[i][j] = prob * 2.0 * volume[i][j] / sum;
				lambda_nf[i][j] = prob * square_side[i] / sum;

			}
		}
	}

void Geometry::Print(int Mx, int My) {

	for (int i = 0; i < Mx; ++i) {
		for (int j = 0; j < My; ++j) {
			cout << lambda_bb[i][j] << " ";
		}
		cout << endl;
	}
	}

	void Geometry::AllocateMemory(int x, int y) {
		lambda_bb.assign(x+2, vector<double>(y+2, 0));
		lambda_bn.assign(x+2, vector<double>(y+2, 0));
		lambda_bf.assign(x+2, vector<double>(y+2, 0));
		lambda_fb.assign(x+2, vector<double>(y+2, 0));
		lambda_fn.assign(x+2, vector<double>(y+2, 0));
		lambda_ff.assign(x+2, vector<double>(y+2, 0));
		lambda_nb.assign(x+2, vector<double>(y+2, 0));
		lambda_nn.assign(x+2, vector<double>(y+2, 0));
		lambda_nf.assign(x+2, vector<double>(y+2, 0));	
		
	};

	
	void Polar::MemoryVectors() {

		volume.assign(Mx + 2, vector<double>(My + 2, 0));;
		square_side.assign(Mx+2,0);
		square_up.assign(Mx+2,0);
		square_front.assign(Mx+2,0);

	}
	
	void Polar::GetValue(int Mx, int My) {
		this->Mx = Mx;
		this->My = My;
	}

	void Polar::UpdateSquareFront() {
		for (int i = 1; i < Mx + 1; ++i)
			square_front[i] = 0.5*My*(i *i - (i - 1) * (i - 1));
	}

	void Polar::UpdateSquareSide() {
		int delta_h = 1;
		for (int i = 1; i < Mx + 1; ++i)
			square_side[i] = delta_h * (i *i - (i - 1) * (i - 1));
	}

	void Polar::UpdateSquareUp() {
		int delta_h = 1;
		for (int i = 1; i < Mx + 1; ++i)
			square_up[i] = delta_h * i * My;
	}

	void Polar::UpdateVolume() {
		int delta_h = 1;
		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j) {
				volume[i][j] = 0.5*My*(i *i - (i - 1) * (i - 1))*delta_h;
			}
		}
	}



	void System::GetValue(int Mx, int My, int N, double theta, double xmin,
		double xmax, double ymin, double ymax) {
		this->Mx = Mx;
		this->My = My;
		this->N = N;
		this->theta = theta;
		this->xmax = xmax;
		this->xmin = xmin;
		this->ymin = ymin;		
		this->ymax = ymax;
	};


	void System::SingularMatrix() {
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < M; ++j) {
				if (i == j)
					A[i][j] = 1;
				else A[i][j] = 0;
			}
		}
	};


	void SCF::MemoryVectors() {
		grad.assign(M+2, 0);
		direction.assign(M + 2, 0);
		alpha.assign(M + 2, 0);
		beta.assign(M+2, 0);
		u.assign(M + 2, 0);
		fi_p.assign(Mx + 2, vector<double>(My + 2, 0));
		fi_w.assign(Mx + 2, vector<double>(My + 2, 0));
		G.assign(Mx + 2, vector<double>(My + 2, 0));
		Gforw.assign(Mx + 1, vector <vector<double>>(My + 1, vector<double>(N + 1, 0)));
		Gback.assign(Mx + 1, vector <vector<double>>(My + 1, vector<double>(N + 1, 0)));
	};

	
	vector <vector <double>> SCF::FindG() {

		for (int i = 0; i < Mx + 1; ++i) {
			for (int j = 0; j < My + 1; ++j) {
				G[i][j] = exp(-u[i*(My + 2) + j]);
			}
		}

		return G;
	};

	vector3d SCF::FindGforw(Geometry& lambda_bb, Geometry& lambda_bn, Geometry& lambda_nb, Geometry& lambda_nn,
		Geometry& lambda_nf, Geometry& lambda_bf, Geometry &lambda_fb, Geometry & lambda_fn, Geometry & lambda_ff) {
		for (int i = xmin; i < xmax + 1; ++i) {
			for (int k = ymin; k < ymax + 1; ++k) {
				Gforw[i][k][0] = G[i][k];
			}
		}

		for (int j = 1; j < N; ++j) {
			for (int k = 1; k < My + 1; ++k) {
				for (int i = 1; i < Mx + 1; ++i) {
					Gforw[i][k][j] = G[i][k] * (lambda_bb[i][k] * Gforw[i - 1][k - 1][j - 1] + lambda_bn[i][k] * Gforw[i - 1][k][j - 1] + lambda_bf[i][k] * Gforw[i - 1][k + 1][j - 1] +
						lambda_nb[i][k] * Gforw[i][k - 1][j - 1] + lambda_nn[i][k] * Gforw[i][k][j - 1] + lambda_nf[i][k] * Gforw[i][k + 1][j - 1] +
						lambda_fb[i][k] * Gforw[i + 1][k - 1][j - 1] + lambda_fn[i][k] * Gforw[i + 1][k + 1][j - 1] + lambda_ff[i][k] * Gforw[i + 1][k + 1][j - 1]);

					Gforw[i][My + 1][j] = Gforw[i][1][j];               
					Gforw[i][0][j] = Gforw[i][My][j];
				}

				Gforw[Mx + 1][k][j] = Gforw[Mx][k][j];                
				Gforw[0][k][j] = 0;
			}
		}


		return Gforw;

	};

	vector3d SCF::FindGback(Geometry& lambda_bb, Geometry& lambda_bn, Geometry& lambda_nb, Geometry& lambda_nn,
		Geometry& lambda_nf, Geometry& lambda_bf, Geometry &lambda_fb, Geometry & lambda_fn, Geometry & lambda_ff) {

		for (int i = 1; i < Mx + 1; ++i) {
			for (int k = 1; k < My + 1; ++k) {
				Gback[i][k][N - 1] = G[i][k];
			}
		}
		for (int j = N - 1; j > 0; --j) {
			for (int k = 1; k < My + 1; ++k) {
				for (int i = 1; i < Mx + 1; ++i) {
					Gback[i][k][j - 1] = G[i][k] * (lambda_bb[i][k] * Gback[i - 1][k - 1][j] + lambda_bn[i][k] * Gback[i - 1][k][j] + lambda_bf[i][k] * Gback[i - 1][k + 1][j] +
						lambda_nb[i][k] * Gback[i][k - 1][j] + lambda_nn[i][k] * Gback[i][k][j] + lambda_nf[i][k] * Gback[i][k + 1][j] +
						lambda_fb[i][k] * Gback[i + 1][k - 1][j] + lambda_fn[i][k] * Gback[i + 1][k][j] + lambda_ff[i][k] * Gback[i + 1][k + 1][j]);

					Gback[i][My + 1][j - 1] = Gback[i][1][j - 1];                
					Gback[i][0][j - 1] = Gback[i][My][j - 1];

				}

				Gback[Mx + 1][k][j - 1] = Gback[Mx + 1][k][j - 1];            
				Gback[0][k][j - 1] = 0;

			}
		}

		return Gback;
	};

	double SCF::FindQ(Polar &volume){

		double sum = 0, q = 0;
		for (int i = 1; i < Mx + 1; ++i) {
			for (int k = 1; k < My + 1; ++k) {
				sum = 0;
				for (int j = 0; j < N; ++j)
					sum += (Gforw[i][k][j] * Gback[i][k][j]) / G[i][k];
				q += sum * volume[i][k];
			}
		}

		return q;
	};

	vector <vector<double>> SCF::FindFiP() {
		double sum, teta = theta * N;
		for (int i = 1; i < Mx + 1; ++i) {
			for (int k = 1; k < My + 1; ++k) {
				sum = 0;
				for (int j = 0; j < N; ++j)
					sum += (Gforw[i][k][j] * Gback[i][k][j]) / G[i][k];

				fi_p[i][k] = (teta / q)*sum;
			}
		}
		return fi_p;
	};

	vector <vector<double>> SCF::FindFiW() {

		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j)
				fi_w[i][j] = G[i][j];
		}

		return fi_w;
	};

	vector <double> SCF::FindGrad() {
		vector <double> grad(M + 1, 0);
		for (int i = 1; i < Mx + 1; ++i) {
			for (int k = 1; k < My + 1; ++k)
				grad[i * (My + 2) + k] = -1 + 1.0 / (fi_p[i][k] + fi_w[i][k]);
		}
		return grad;
	};
