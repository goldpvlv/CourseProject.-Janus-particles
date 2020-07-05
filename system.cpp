#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include <conio.h>


#include "molecule.h"
#include "method.h"
#include "geometry.h"
#include "system.h"



using namespace std;


System::System(string _source) {
	source = _source;
}

void System::ReadParameters() {
	
	ifstream ifs(source, ifstream::in);
	string word;
	while (ifs.good()) {
		ifs >> word;
		if (word.find("geometry") != string::npos) ifs >> geometry_name;

		if (word.find("layers_x") != string::npos) ifs >> layers_x;
		if (word.find("layers_y") != string::npos) ifs >> layers_y;
		if (word.find("curvature") != string::npos) ifs >> Rcurv;
	}
	ifs.close();
}

void System::SetGeometry() {

	if (geometry_name == "polar") {
		Polar polar;
		polar.GetValue(layers_x, layers_y, Rcurv);
		polar.AllocateMemory();
		polar.UpdateVolume();
		polar.UpdateSquareFront();
		polar.UpdateSquareUp();
		polar.UpdateSquareSide();
		polar.Transposition();
		
		geo = polar;
	}
	else if (geometry_name == "torus") {
		Torus torus;
		torus.GetValue(layers_x, layers_y, Rcurv);
		torus.AllocateMemory();
		torus.UpdateVolume();
		torus.UpdateSquareFront();
		torus.UpdateSquareUp();
		torus.UpdateSquareLeft();
		torus.UpdateSquareRight();
		torus.Transposition();

		geo = torus;
	}
	else {
		printf("ERROR: unknown geometry\n");
		return;
	}
	
	//geo.PrintVolume();

}

void System::ReadMolecules() {
	
	ifstream ifs(source, ifstream::in);
	string word;
	while (ifs.good()) {
		ifs >> word;
		if (word.find("molecule") != string::npos) {
			Molecule new_mol;
			while (word != "[" && ifs.good()) {
				ifs >> word;
				if (word.find("Ns") != string::npos) ifs >> new_mol.ns;
				if (word.find("gen") != string::npos) ifs >> new_mol.num_generation;
				if (word.find("sigma") != string::npos) ifs >> new_mol.sigma;
				if (word.find("chi") != string::npos) ifs >> new_mol.chi;
				if (word.find("xmin") != string::npos) ifs >> new_mol.xmin;
				if (word.find("xmax") != string::npos) ifs >> new_mol.xmax;
				if (word.find("ymin") != string::npos) ifs >> new_mol.ymin;
				if (word.find("ymax") != string::npos) ifs >> new_mol.ymax;

			}

			new_mol.MyNewMethod();
			mol.push_back(new_mol);
		}
	}
	ifs.close();

	
	for (int i = 0; i < mol.size(); i++) {
		cout << "number of molecule: " << i + 1 << endl;
		cout << mol[i].ns << endl;
		cout << mol[i].num_generation << endl;
		cout << mol[i].xmin << endl;

	}
};

void System::ReadMethods() {

	ifstream ifs(source, ifstream::in);
	string word;
	while (ifs.good()) {
		ifs >> word;
		if (word.find("method") != string::npos) {
			double tolerance_, nu_, _num_iter;
			string name_method;
			// пока не встретим новый блок или не закончитс€ файл
			while (word != "[" && ifs.good()) {
				ifs >> word;
				if (word.find("type") != string::npos) ifs >> name_method;
				if (word.find("tolerance") != string::npos) ifs >> tolerance_;
				if (word.find("num_iter") != string::npos) ifs >> _num_iter;
				if (word.find("step") != string::npos) ifs >> nu_;
			}
			//num_iter_, tolerance_, step_
			if (name_method == "gradient") {
				methods.push_back(make_unique<Gradient>(_num_iter, tolerance_, nu_));
			}
			else if (name_method == "DFP") {
				methods.push_back(make_unique<DFP>(_num_iter, tolerance_, nu_));
			}
			else {
				printf("ERROR: unknown method");
			}

		}
	}
	ifs.close();

	/*for (BaseOptimTools &w : methods) {
		cout << w.name << endl;
		cout << w.tolerance << endl;
		cout << w.nu << endl;
		cout << w.num_iter << endl;
	}*/
}

void System::Function() {

	int M = 3 * (layers_x + 2)*(layers_y + 2);

	vector < vector < double > > fi_solv(layers_x + 1, vector<double>(layers_y + 1, 0));
	vector < vector < double > > fi_solv(layers_x + 1, vector<double>(layers_y + 1, 0));
	vector <double> u_solv(M + 2, 0);


	for (int t = 0; t < mol.size(); ++t) {

		mol[t].fi_p.assign(layers_x + 2, vector<double>(layers_y + 2, 0));
		mol[t].fi_w.assign(layers_x + 2, vector<double>(layers_y + 2, 0));
		mol[t].G.assign(layers_x + 2, vector<double>(layers_y + 2, 0));
		mol[t].GSide.assign(M + 2);
		mol[t].Gforw.assign(layers_x + 2, vector <vector<double>>(layers_y + 2, vector<double>(mol[t].ns + 2, 0)));
		mol[t].Gback.assign(layers_x + 2, vector <vector<double>>(layers_y + 2, vector<double>(mol[t].ns + 2, 0)));

		FindG(t);
		FindGSide(t);

		for (int i = mol[t].xmin; i < mol[t].xmax + 1; ++i) {
			for (int k = mol[t].ymin; k < mol[t].ymax + 1; ++k) {
				mol[t].Gforw[i][k][0] = mol[t].G[i][k];
			}
		}
		FindGforw(t);

		for (int i = 1; i < layers_x + 1; ++i) {
			for (int k = 1; k < layers_y + 1; ++k) {
				mol[t].Gback[i][k][mol[t].ns - 1] = mol[t].G[i][k];
			}
		}
		FindGback(t);

		FindQ(t);
		FindFiP(t);
		FindFiW(t);
		FindLagrangeMultipliers(t);
	}

	FindGSolv();
	FindGTotal();
	FindLagrangeMultipliers();
		
		grad = FindGrad();

	

};

void System::FindG(int t) {


	for (int i = 0; i < layers_x + 1; ++i) {
		for (int j = 0; j < layers_y + 1; ++j) {
			mol[t].G[i][j] = exp(-mol[t].u[i*(layers_y + 2) + j]);
		}
	}

};

void System::FindGSide(int t) {

	double sum = 0.0;

		for (int i = -1; i <= 1; ++i) {
			for (int j = -1; j <= 1; ++j) {
				if (i == j == -1) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_bb[l][k] * mol[t].G[l][k];
					}
				}
				else if (i == -1 && j == 0) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_bn[l][k] * mol[t].G[l][k];
					}
				}
				else if (i == -1 && j == 1) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_bf[l][k] * mol[t].G[l][k];
					}
				}
				else if (i == 0 && j == -1) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_nb[l][k] * mol[t].G[l][k];
					}
				}
				else if (i == 0 && j == 0) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_nn[l][k] * mol[t].G[l][k];
					}
				}
				else if (i == 0 && j == 1) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_nf[l][k] * mol[t].G[l][k];
					}
				}
				else if (i == 1 && j == -1) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_fb[l][k] * mol[t].G[l][k];
					}
				}
				else if (i == 1 && j == 0) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_fn[l][k] * mol[t].G[l][k];
					}
				}
				else if (i == 1 && j == 1) {
					for (int l = 0; l < layers_x + 1; ++l) {
						for (int k = 0; k < layers_y + 1; ++k)
							mol[t].GSide[l * (layers_y + 2) + k] = geo.lambda_ff[l][k] * mol[t].G[l][k];
					}
				}
			}

		}
}


void System::FindGforw(int t) {

	vector3d Gforw(layers_x + 2, vector <vector<double>>(layers_y + 2, vector<double>(mol[t].ns, 0)));

	for (int i = mol[t].xmin; i < mol[t].xmax + 1; ++i) {
		for (int k = mol[t].ymin; k < mol[t].ymax + 1; ++k) {
			mol[t].Gforw[i][k][0] = mol[t].G[i][k];
		}
	}

	for (int j = 1; j < mol[t].ns; ++j) {
		for (int k = 1; k < layers_x + 1; ++k) {
			for (int i = 1; i < layers_y + 1; ++i) {
				mol[t].Gforw[i][k][j] = mol[t].G[i][k] * (geo.lambda_bb[i][k] * mol[t].Gforw[i - 1][k - 1][j - 1] + geo.lambda_bn[i][k] * mol[t].Gforw[i - 1][k][j - 1] + geo.lambda_bf[i][k] * mol[t].Gforw[i - 1][k + 1][j - 1] +
					geo.lambda_nb[i][k] * mol[t].Gforw[i][k - 1][j - 1] + geo.lambda_nn[i][k] * mol[t].Gforw[i][k][j - 1] + geo.lambda_nf[i][k] * mol[t].Gforw[i][k + 1][j - 1] +
					geo.lambda_fb[i][k] * mol[t].Gforw[i + 1][k - 1][j - 1] + geo.lambda_fn[i][k] * mol[t].Gforw[i + 1][k + 1][j - 1] + geo.lambda_ff[i][k] * mol[t].Gforw[i + 1][k + 1][j - 1]);

				mol[t].Gforw[i][layers_y + 1][j] = mol[t].Gforw[i][1][j];
				mol[t].Gforw[i][0][j] = mol[t].Gforw[i][layers_x][j];
			}

			mol[t].Gforw[layers_x + 1][k][j] = mol[t].Gforw[layers_x][k][j];
			mol[t].Gforw[0][k][j] = 0;
		}
	}

};


void System::FindGback(int t) {

	for (int i = 1; i < layers_x + 1; ++i) {
		for (int k = 1; k < layers_y + 1; ++k) {
			mol[t].Gback[i][k][mol[t].ns - 1] = mol[t].G[i][k];
		}
	}
	for (int j = mol[t].ns - 1; j > 0; --j) {
		for (int k = 1; k < layers_y + 1; ++k) {
			for (int i = 1; i < layers_x + 1; ++i) {
				mol[t].Gback[i][k][j - 1] = mol[t].G[i][k] * (geo.lambda_bb[i][k] * mol[t].Gback[i - 1][k - 1][j] + geo.lambda_bn[i][k] * mol[t].Gback[i - 1][k][j] + geo.lambda_bf[i][k] * mol[t].Gback[i - 1][k + 1][j] +
					geo.lambda_nb[i][k] * mol[t].Gback[i][k - 1][j] + geo.lambda_nn[i][k] * mol[t].Gback[i][k][j] + geo.lambda_nf[i][k] * mol[t].Gback[i][k + 1][j] +
					geo.lambda_fb[i][k] * mol[t].Gback[i + 1][k - 1][j] + geo.lambda_fn[i][k] * mol[t].Gback[i + 1][k][j] + geo.lambda_ff[i][k] * mol[t].Gback[i + 1][k + 1][j]);

				mol[t].Gback[i][layers_y + 1][j - 1] = mol[t].Gback[i][1][j - 1];
				mol[t].Gback[i][0][j - 1] = mol[t].Gback[i][layers_y][j - 1];

			}

			mol[t].Gback[layers_x + 1][k][j - 1] = mol[t].Gback[layers_x + 1][k][j - 1];
			mol[t].Gback[0][k][j - 1] = 0;

		}
	}
};


void System::FindQ(int t) {

	double sum = 0, q = 0;
	for (int i = 1; i < layers_x + 1; ++i) {
		for (int k = 1; k < layers_y + 1; ++k) {
			sum = 0;
			for (int j = 0; j < mol[t].ns; ++j)
				sum += (mol[t].Gforw[i][k][j] * mol[t].Gback[i][k][j]) / mol[t].G[i][k];
			q += sum * geo.volume[i][k];
		}
	}

};

void System::FindFiP(int t) {
	double sum;
	for (int i = 1; i < layers_x + 1; ++i) {
		for (int k = 1; k < layers_y + 1; ++k) {
			sum = 0;
			for (int j = 0; j < mol[t].ns; ++j)
				sum += (mol[i].Gforw[i][k][j] * mol[i].Gback[i][k][j]) / G[i][k];

			mol[t].fi_p[i][k] = (mol[t].theta / mol[t].q)*sum;
		}
	}
};

void System::FindFiW(int t) {

	for (int i = 1; i < layers_x + 1; ++i) {
		for (int j = 1; j < layers_y + 1; ++j)
			mol[t].fi_w[i][j] = mol[t].G[i][j];
	}

};

void System::FindLagrangeMultipliers(int t) {

	mol[t].u

};

void System::FindGSolv() {

	for (int i = 1; i < layers_x + 1; ++i) {
		for (int j = 1; j < layers_y + 1; ++j)
			G_solv[i][j]= exp(-u_solv[i*(layers_y + 2) + j]);
	}

}

void System::FindGTotal() {

	for (int i = 1; i < layers_x + 1; ++i) {
		for (int j = 1; j < layers_y + 1; ++j)
			G_total[i][j] = mol[1].G[i][j] + mol[2].G[i][j] + G_solv[i][j];
	}
}

void System::Cycling() {
	for (BaseOptimTools &scheme : methods) {
		cout << scheme.name << endl;


		scheme.SetParameters(layers_x, layers_y, grad, u);

		double length_of_grad = 0.0;
		int step = 0.0;

		//Function();
		length_of_grad = scheme.SetGradFirst(grad);

		do {
			scheme.UpdateX(u, grad);
			//Function();
			length_of_grad = scheme.SetGradRegular(grad);

			step++;
			cout << step << "     " << length_of_grad << endl;
		} while ((length_of_grad > scheme.tolerance) && (step < scheme.num_iter));
	}
};