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

	void Polar::PrintVolume() {
		for (int i = 1; i < Mx + 1; ++i)
			cout << square_front[i] << " ";
	};


	void Polar::TranspositionPolar(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn,
		vector<vector<double>> & lambda_nb, vector<vector<double>> & lambda_nn, vector<vector<double>> & lambda_nf,
		vector<vector<double>> & lambda_bf, vector<vector<double>> &lambda_fb, vector<vector<double>> & lambda_fn, vector<vector<double>> & lambda_ff) {
		for (int i = 1; i < Mx + 1; ++i) {
			for (int j = 1; j < My + 1; ++j) {
				lambda_bb[i][j] = 0.0;
				lambda_bn[i][j] = square_front[i - 1] / volume[i][j] * 1.0 / 6.0;
				cout << square_front[i - 1] << " "<<volume[i][j] << " " << lambda_bb[i][j]<<endl;
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



	void System::GetValue(int Mx, int My, int N, int M, double theta, double xmin,
		double xmax, double ymin, double ymax) {
		this->M = M;
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
		A.assign(M, vector<double>(M, 0));
		grad.assign(M+2, 0);
		direction.assign(M + 2, 0);
		alpha.assign(M + 2, 0);
		beta.assign(M+2, 0);
		u.assign(M + 2, 0);
		fi_p.assign(Mx + 2, vector<double>(My + 2, 0));
		fi_w.assign(Mx + 2, vector<double>(My + 2, 0));
		G.assign(Mx + 2, vector<double>(My + 2, 0));
		Gforw.assign(Mx + 2, vector <vector<double>>(My + 2, vector<double>(N + 2, 0)));
		Gback.assign(Mx + 2, vector <vector<double>>(My + 2, vector<double>(N + 2, 0)));
	};

	
	vector <vector <double>> SCF::FindG() {

		for (int i = 0; i < Mx + 1; ++i) {
			for (int j = 0; j < My + 1; ++j) {
				G[i][j] = exp(-u[i*(My + 2) + j]);
			}
		}

		return G;
	};

	vector3d SCF::FindGforw(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn, vector<vector<double>> & lambda_nb, 
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

	vector3d SCF::FindGback(vector<vector<double>> & lambda_bb, vector<vector<double>> & lambda_bn, vector<vector<double>> & lambda_nb, 
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

	double SCF::FindQ(vector<vector<double>> & volume){

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


	void SCF::Print(){
		for (int i = 1; i < M + 1; ++i) {
			cout << grad[i] << " ";
		}

	};
	vector < double >  product_matrix_on_vector(vector < vector < double > >A, vector <double> grad, int M) {
		vector <double> a(M, 0);

		for (int i = 0; i < M; i++) {
			a[i] = 0;
			for (int j = 0; j < M; j++) {
				a[i] += A[i][j] * grad[j];
			}
		}

		return a;
	}

	vector < double >  System::FindDirection() {
	
		direction = product_matrix_on_vector(A, grad, M);

		for (int i = 0; i < M; i++)
			direction[i] = -1 * direction[i];
		return direction;
	
	};

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


