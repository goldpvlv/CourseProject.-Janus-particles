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


void Geometry::Print() {

	for (int i = 0; i < Mx; ++i) {
		for (int j = 0; j < My; ++j) {
			cout << lambda_bb[i][j] << " ";
		}
		cout << endl;
	}
	}

	void Geometry::AllocateMemory() {
		lambda_bb.assign(Mx+2, vector<double>(My+2, 0));
		lambda_bn.assign(Mx+2, vector<double>(My+2, 0));
		lambda_bf.assign(Mx+2, vector<double>(My+2, 0));
		lambda_fb.assign(Mx+2, vector<double>(My+2, 0));
		lambda_fn.assign(Mx+2, vector<double>(My+2, 0));
		lambda_ff.assign(Mx+2, vector<double>(My+2, 0));
		lambda_nb.assign(Mx+2, vector<double>(My+2, 0));
		lambda_nn.assign(Mx+2, vector<double>(My+2, 0));
		lambda_nf.assign(Mx+2, vector<double>(My+2, 0));	
		volume.assign(Mx + 2, vector<double>(My + 2, 0));;
		square_side.assign(Mx + 2, vector<double>(My + 2, 0));
		square_up.assign(Mx + 2, vector<double>(My + 2, 0));
		square_front.assign(Mx + 2, vector<double>(My + 2, 0));
		square_right.assign(Mx + 2, vector<double>(My + 2, 0));;
		square_left.assign(Mx + 2, vector<double>(My + 2, 0));;
		
	};

	
	void Geometry::GetValue(int Mx, int My, int N) {
		this->Mx = Mx;
		this->My = My;
		this->N = N;
	}


	void Polar::PrintVolume() {
		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j) {
				cout << volume[i][j] << " ";
			}
			cout << endl;
		}
	};


	void Polar::Transposition(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn,
		vector<vector<double>> & lambda_nb, vector<vector<double>> & lambda_nn, vector<vector<double>> & lambda_nf,
		vector<vector<double>> & lambda_bf, vector<vector<double>> &lambda_fb, vector<vector<double>> & lambda_fn,
		vector<vector<double>> & lambda_ff, vector<vector<double>> & volume, vector<vector<double>> & square_up,
		vector<vector<double>> & square_side, vector<vector<double>> & square_front) {

		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j) {
				lambda_bb[i][j] = 0.0;
				lambda_bn[i][j] = square_front[i - 1][j] / volume[i][j] * 1.0 / 6.0;
				cout << square_front[i - 1][j] << " " << volume[i][j] << " " << lambda_bb[i][j] << endl;
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
		int delta_h = 1;
		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j) {
				square_side[i][j] = delta_h * (i *i - (i - 1) * (i - 1));
			}
		}
	}

	void Polar::UpdateSquareUp() {
		int delta_h = 1;
		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j) {
				square_up[i][j] = delta_h * i * My;
			}
		}
	}

	void Polar::UpdateVolume() {
		int delta_h = 1;
		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j) {
				volume[i][j] = 0.5*My*(i *i - (i - 1) * (i - 1))*delta_h;
			}
		}
	}

	void Torus::UpdateSquareFront() {


		for (int i = 1; i < Mx + 1; i++)
			for (int j = 1; j < My + 1; j++) {
				square_front[i][j] = i * delta_alpha * 2.0 * pi*(N + i * cos((2 * j - 1)*delta_alpha) / 2.0);
			}

	}

	void Torus::UpdateVolume() {
		vector <double> rho(Mx + 1, 0);

		for (int i = 1; i < Mx + 1; ++i) {
			rho[i] = 0.67*(2 * i - 1);
		}

		for (int i = 1; i < Mx + 1; i++) {
			for (int j = 1; j < My + 1; j++) {
				volume[i][j] = square_up[i][j] * pi * 2.0*(N + rho[i] * cos(((2 * j - 1) / 2.0)*delta_alpha));
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
				square_right[i][j] = (i - (i - 1)) * 2.0 * pi*(N + ((2 * i - 1) / 2.0 * cos(j*delta_alpha)));
			}

	};

	void Torus::UpdateSquareLeft() {

		for (int i = 1; i < Mx + 1; i++)
			for (int j = 1; j < My + 1; j++) {
				square_left[i][j] = (i - (i - 1)) * 2.0 * pi*(N + ((2 * i - 1) / 2.0*cos((j - 1)*delta_alpha)));
			}

	};

	void Torus::Transposition(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn,
		vector<vector<double>> & lambda_nb, vector<vector<double>> & lambda_nn, vector<vector<double>> & lambda_nf,
		vector<vector<double>> & lambda_bf, vector<vector<double>> &lambda_fb, vector<vector<double>> & lambda_fn,
		vector<vector<double>> & lambda_ff, vector<vector<double>> & volume, vector<vector<double>> & square_up,
		vector<vector<double>> & square_left, vector<vector<double>> & square_right,  vector<vector<double>> & square_front) {
	
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


	void Molecule::GetValue(int Mx, int My, int N, double theta, double xmin, double xmax, double ymin, double ymax) {

		this->Mx = Mx;
		this->My = My;
		this->N = N;
		this->theta = theta;
		this->xmax = xmax;
		this->xmin = xmin;
		this->ymin = ymin;		
		this->ymax = ymax;
	};


	void Molecule::MemoryVectors() {

		fi_p.assign(Mx + 2, vector<double>(My + 2, 0));
		fi_w.assign(Mx + 2, vector<double>(My + 2, 0));
		G.assign(Mx + 2, vector<double>(My + 2, 0));
		Gforw.assign(Mx + 2, vector <vector<double>>(My + 2, vector<double>(N + 2, 0)));
		Gback.assign(Mx + 2, vector <vector<double>>(My + 2, vector<double>(N + 2, 0)));
	};

	
	vector <vector <double>> Molecule::FindG(vector <double> & u) {

		for (int i = 0; i < Mx + 1; ++i) {
			for (int j = 0; j < My + 1; ++j) {
				G[i][j] = exp(-u[i*(My + 2) + j]);
			}
		}

		return G;
	};

	vector3d Molecule::FindGforw(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn, vector<vector<double>> & lambda_nb, 
		vector<vector<double>> & lambda_nn, vector<vector<double>> & lambda_nf, vector<vector<double>> & lambda_bf, 
		vector<vector<double>> &lambda_fb, vector<vector<double>> & lambda_fn, vector<vector<double>> & lambda_ff) {
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

	vector3d Molecule::FindGback(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn, vector<vector<double>> & lambda_nb, 
		vector<vector<double>> & lambda_nn, vector<vector<double>> & lambda_nf, vector<vector<double>> & lambda_bf, 
		vector<vector<double>> & lambda_fb, vector<vector<double>> & lambda_fn, vector<vector<double>> & lambda_ff) {

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

	double Molecule::FindQ(vector<vector<double>> & volume){

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

	vector <vector<double>> Molecule::FindFiP() {
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

	vector <vector<double>> Molecule::FindFiW() {

		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j)
				fi_w[i][j] = G[i][j];
		}

		return fi_w;
	};

	vector <double> System::FindGrad(vector < vector < double > > &fi_p, vector < vector < double > > &fi_w) {
		vector <double> grad(M + 1, 0);
		for (int i = 1; i < Mx + 1; ++i) {
			for (int k = 1; k < My + 1; ++k)
				grad[i * (My + 2) + k] = -1 + 1.0 / (fi_p[i][k] + fi_w[i][k]);
		}
		return grad;
	};


	

	void System::GetValue(int Mx, int My) {
		this -> Mx = Mx;
		this -> My = My;
	
	};


	void System::MemoryVectors() {
		grad.assign(M + 2, 0);
		u.assign(M + 2, 0);
		direction.assign(M + 2, 0);
		alpha.assign(M + 2, 0);
		beta.assign(M + 2, 0);
		A.assign(M, vector<double>(M, 0));
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

	vector < double >  System::FindDirection() {
	
		direction = product_matrix_on_vector(A, grad, M);

		for (int i = 0; i < M; i++)
			direction[i] = -1 * direction[i];
		return direction;
	
	};


	vector < vector < double > > System::Formula() {
		double number1 = 0, number2 = 0;
		vector < vector < double > > matrix1(M, vector<double>(M, 0));
		vector < vector < double > > matrix2(M, vector<double>(M, 0));
		vector < vector < double > > B(M, vector<double>(M, 0));
		vector < vector < double > > C(M, vector<double>(M, 0));
		vector <double> a(M, 0);
		vector <double> b(M, 0);
		matrix1 = make_matrix(alpha, alpha, M);
		number1 = make_number(alpha, beta, M);
		B = division_matrix_by_number(matrix1, number1, M);
		a = product_matrix_on_vector(A, beta, M);
		matrix2 = make_matrix(a, beta, M);
		matrix2 = product(matrix2, A, M);
		b = product_matrix_on_vector(A, beta, M);
		number2 = make_number(b, beta, M);
		C = division_matrix_by_number(matrix2, number2, M);

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++)
				A[i][j] = A[i][j] + B[i][j] - C[i][j];
		}

		return A;

	};

	void System::LengthOfGrad() {

		double sum = 0;
		for (int i = 0; i < M; ++i) {
			sum += grad[i] * grad[i];
		}
		length_of_grad = sqrt(sum);

	};


	void System::Print(){
			for (int i = 1; i < M + 1; ++i) {
				cout << grad[i] << " ";
			}

		};


	vector < double >  product_matrix_on_vector(vector < vector < double > >A, vector <double> q, int M) {
		vector <double> a(M, 0);

		for (int i = 0; i < M; i++) {
			a[i] = 0;
			for (int j = 0; j < M; j++) {
				a[i] += A[i][j] * q[j];
			}
		}

		return a;
	}

	vector < vector < double > > product(vector < vector < double > > A, vector < vector < double > > U, int N) {
		vector < vector < double > > c(N, vector<double>(N, 0));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				c[i][j] = 0;
				for (int t = 0; t < N; t++)
					c[i][j] += A[i][t] * U[t][j];
			}
		}
		return c;
	}

	vector < vector < double > > make_matrix(vector <double> a, vector <double> b, int M) {
		vector < vector < double > > matrix(M, vector<double>(M, 0));
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++)
				matrix[i][j] = a[i] * b[j];
		}
		return matrix;
	}


	double make_number(vector <double> a, vector <double> b, int N) {
		double number = 0;
		for (int i = 0; i < N; i++)
			number += a[i] * b[i];
		return number;
	}



	vector < vector < double > > division_matrix_by_number(vector < vector < double > > matrix, double number, int N) {
		vector < vector < double > > matrix1(N, vector<double>(N, 0));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++)
				matrix1[i][j] = matrix[i][j] / number;
		}
		return matrix1;
	}