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
const double pi = 3.14159265359;


void DFP::SingularMatrix() {
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < M; ++j) {
			if (i == j)
				A[i][j] = 1;
			else A[i][j] = 0;
		}
	}
};




vector < vector < double > > DFP::Formula() {
	double number1 = 0, number2 = 0;
	vector < vector < double > > matrix1(M, vector<double>(M, 0));
	vector < vector < double > > matrix2(M, vector<double>(M, 0));
	vector < vector < double > > B(M, vector<double>(M, 0));
	vector < vector < double > > C(M, vector<double>(M, 0));

	vector <double> a(M, 0);
	vector <double> b(M, 0);


	matrix1 = Multiply_vectors(alpha, alpha);
	number1 = Find_number(alpha, beta);
	B = Division_matrix_on_number(matrix1,number1);
	a = Multiply_matrix_by_vector(A,beta);
	matrix2 = Multiply_vectors(a, beta);
	matrix2 = Multiply_matrixes(A, matrix1);
	b = Multiply_matrix_by_vector(A, beta);
	number2 = Find_number(b, beta);
	C = Division_matrix_on_number(matrix2 ,number2);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++)
			A[i][j] = A[i][j] + B[i][j] - C[i][j];
	}

	return A;

};

vector < vector < double > > DFP::Multiply_matrixes(vector < vector < double > > A, vector < vector < double > > U) {

	vector < vector < double > > c(M, vector<double>(M, 0));

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			c[i][j] = 0;
			for (int t = 0; t < M; t++)
				c[i][j] += A[i][t] * U[t][j];
		}
	}
	return c;
};


vector < vector < double > > DFP::Multiply_vectors(vector<double> a, vector<double> b) {

	vector < vector < double > > matrix(M, vector<double>(M, 0));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++)
			matrix[i][j] = a[i] * b[j];
	}
	return matrix;
};

double DFP::Find_number(vector<double> a , vector<double> b ) {
	double number = 0;
	for (int i = 0; i < M; i++)
		number += a[i] * b[i];
	return number;
};


 vector < double > DFP::Multiply_matrix_by_vector(vector<vector<double>>A,vector<double>grad ) {
	vector <double> a(M, 0);
	for (int i = 0; i < M; i++) {
		a[i] = 0;
		for (int j = 0; j < M; j++) {
			a[i] += A[i][j] * grad[j];
		}
	}
	return a;
};

 vector < double > BaseOptimTools::Multiply_matrix_by_vector(vector<vector<double>>A, vector<double>grad) {
	 vector <double> a(M, 0);
	 for (int i = 0; i < M; i++) {
		 a[i] = 0;
		 for (int j = 0; j < M; j++) {
			 a[i] += A[i][j] * grad[j];
		 }
	 }
	 return a;
 };

 vector < vector < double > > DFP::Division_matrix_on_number(vector<vector<double>>A, double b) {
	 vector < vector < double > > matrix1(M, vector<double>(M, 0));
	 for (int i = 0; i < M; i++) {
		 for (int j = 0; j < M; j++)
			 matrix1[i][j] = A[i][j] / b;
	 }
	 return matrix1;
 };

//BaseOptimTools operator/(int a, const BaseOptimTools & matrix)  //division_matrix_by_number
//{
//	BaseOptimTools tmp(int M);
//	for (int i = 0; i < M; i++)
//		for (int j = 0; j < M; j++)
//
//			tmp[i][j] += a * matrix[i][j];
//	return tmp;
//}
//
//BaseOptimTools operator*(const BaseOptimTools & matrix) {// multiply 2 matrixes
//
//	BaseOptimTools temp(matrix.M);
//	temp.M = matrix.M;
//	 for (int i = 0; i < this->M; i++)
//	 {
//		 double tm = 0;
//		 for (int j = 0; j < temp.M; j++)
//		 {
//			 for (int k = 0; k < temp.M; k++)
//				 tm += this->matrix[i][k] * matrix.mat[k][j];
//			 temp.matrix[i][j] = tm;
//			 tm = 0;
//		 }
//	 }
// return temp;
//};  



double Gradient::LengthOfGrad() {

	double sum = 0;
	for (int i = 0; i < M; ++i) {
		sum += grad[i] * grad[i];
	}
	length_of_grad = sqrt(sum);

	return length_of_grad;

};


void BaseOptimTools::Print() {
	for (int i = 1; i < M + 1; ++i) {
		cout << grad[i] << " ";
	}

};

vector<double> BaseOptimTools::FindDirection(vector<vector<double>>A, vector<double>grad) {
	
	vector <double> direction(M, 0);

	direction = Multiply_matrix_by_vector(A, grad);

	for (int i = 0; i < M; i++)
		direction[i] = -1 * direction[i];
	return direction;

};

vector <double> BaseOptimTools::Function(){

	vector3d Gback(Mx + 2, vector <vector<double>>(My + 2, vector<double>(R0, 0)));
	vector3d Gforw(Mx + 2, vector <vector<double>>(My + 2, vector<double>(R0, 0)));
	vector <vector <double>> G(Mx + 2, vector <double>(My + 2, 0));
	vector <double> grad(M, 0);



	G = FindG();

	for (int i = xmin; i < xmax + 1; ++i) {
		for (int k = ymin; k < ymax + 1; ++k) {
			Gforw[i][k][0] = G[i][k];
		}
	}
	Gforw = FindGforw();


	for (int i = 1; i < Mx + 1; ++i) {
		for (int k = 1; k < My + 1; ++k) {
			Gback[i][k][R0 - 1] = G[i][k];
		}
	}
	Gback = FindGforw();

	q = FindQ();;
	fi_p = FindFiP();
	fi_w = FindFiW();
	grad = FindGrad();

	return grad;

};


void BaseOptimTools::Сycling() {

	SingularMatrix();
	length_of_grad = LengthOfGrad();
	grad = Function();


	do {
		direction = FindDirection(A, grad);

		for (int i = 0; i < M; i++)
			alpha[i] = u[i];

		for (int i = 0; i < M; i++)
			u[i] = u[i] + nu * direction[i];

		for (int i = 0; i < M; i++)
			alpha[i] = u[i] - alpha[i];

		for (int i = 0; i < M; i++)
			beta[i] = grad[i];

		grad = Function();

		for (int i = 0; i < M; i++)
			beta[i] = grad[i] - beta[i];

		A = Formula();                          //Ak=Ak-1 + B - C;

		step++;

		cout << step << "     " << length_of_grad << endl;
	} while ((length_of_grad > tolerance) && (step < num_iter));


};
Gradient::Gradient(double t) {
	name = "Gradient"; 
	tolerance = t;
}

DFP::DFP(double t) {
	name = "DFP"; 
	tolerance = t;
}

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
}


void Geometry::PrintVolume() {
	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j) {
			cout << volume[i][j] << " ";
		}
		cout << endl;
	}
};


void Polar::Transposition() {

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


BaseOptimTools::BaseOptimTools(int M) {
	arr.assign(M + 2, vector<double>(M + 2, 0));
};

void BaseOptimTools::GetValue(int Mx_, int My_, int R0_) {
	Mx = Mx_;
	My = My_;
	R0 = R0_;
}


void BaseOptimTools::AllocateMemory() {

	fi_p.assign(Mx + 2, vector<double>(My + 2, 0));
	fi_w.assign(Mx + 2, vector<double>(My + 2, 0));
	G.assign(Mx + 2, vector<double>(My + 2, 0));
	Gforw.assign(Mx + 2, vector <vector<double>>(My + 2, vector<double>(R0 + 2, 0)));
	Gback.assign(Mx + 2, vector <vector<double>>(My + 2, vector<double>(R0 + 2, 0)));
	A.assign(M + 2, vector<double>(M + 2, 0));
	p. assign(M+2,0);
	grad. assign(M+2,0);
	direction. assign (M+2,0);
	alpha.assign(M+2,0);
	beta.assign(M+2,0);
	u.assign(M+2,0);
};


vector <vector <double>> Gradient::FindG() {

	for (int i = 0; i < Mx + 1; ++i) {
		for (int j = 0; j < My + 1; ++j) {
			G[i][j] = exp(-u[i*(My + 2) + j]);
		}
	}

	return G;
};

vector3d Gradient::FindGforw() {
	for (int i = xmin; i < xmax + 1; ++i) {
		for (int k = ymin; k < ymax + 1; ++k) {
			Gforw[i][k][0] = G[i][k];
		}
	}

	for (int j = 1; j < R0; ++j) {
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

vector3d Gradient::FindGback() {

	for (int i = 1; i < Mx + 1; ++i) {
		for (int k = 1; k < My + 1; ++k) {
			Gback[i][k][R0 - 1] = G[i][k];
		}
	}
	for (int j = R0 - 1; j > 0; --j) {
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

double Gradient::FindQ() {

	double sum = 0, q = 0;
	for (int i = 1; i < Mx + 1; ++i) {
		for (int k = 1; k < My + 1; ++k) {
			sum = 0;
			for (int j = 0; j < R0; ++j)
				sum += (Gforw[i][k][j] * Gback[i][k][j]) / G[i][k];
			q += sum * volume[i][k];
		}
	}

	return q;
};

vector <vector<double>> Gradient::FindFiP() {
	double sum;
	for (int i = 1; i < Mx + 1; ++i) {
		for (int k = 1; k < My + 1; ++k) {
			sum = 0;
			for (int j = 0; j < R0; ++j)
				sum += (Gforw[i][k][j] * Gback[i][k][j]) / G[i][k];

			fi_p[i][k] = (teta / q)*sum;
		}
	}
	return fi_p;
};

vector <vector<double>> Gradient::FindFiW() {

	for (int i = 1; i < Mx + 1; ++i) {
		for (int j = 1; j < My + 1; ++j)
			fi_w[i][j] = G[i][j];
	}

	return fi_w;
};

vector <double> Gradient::FindGrad() {

	vector <double> grad(M + 1, 0);
	for (int i = 1; i < Mx + 1; ++i) {
		for (int k = 1; k < My + 1; ++k)
			grad[i * (My + 2) + k] = -1 + 1.0 / (fi_p[i][k] + fi_w[i][k]);
	}
	return grad;
};

System::System(string _source) {
	source = _source;
}

void System::ReadParameters() {
	// читаем входной файл - ищем ключевое слово - узнаем тип Геометрии
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
	// устанавливаем Геометрию
	if (geometry_name == "polar") {
		Polar polar;
		polar.GetValue(layers_x, layers_y, Rcurv);
		polar.AllocateMemory();
		polar.UpdateVolume();
		// передаем значение полей (вероятности перехода и т.п.) в объект абстрактного класса
		geo = polar;
	}
	else if (geometry_name == "torus") {
		Torus torus;
		torus.GetValue(layers_x, layers_y, Rcurv);
		torus.AllocateMemory();
		torus.UpdateVolume();

		// передаем значение полей (вероятности перехода и т.п.) в объект абстрактного класса
		geo = torus;
	}
	else {
		printf("ERROR: unknown geometry\n");
		return;
	}
	// для теста:
	geo.PrintVolume();

}

void System::ReadMolecules() {
	// читаем входной файл - ищем ключевое слово - узнаем тип Геометрии
	ifstream ifs(source, ifstream::in);
	string word;
	while (ifs.good()) {
		ifs >> word;
		if (word.find("molecule") != string::npos) {
			Molecule new_mol;
			// пока не встретим новый блок или не закончится файл
			while (word != "[" && ifs.good()) {
				ifs >> word;
				if (word.find("Ns") != string::npos) ifs >> new_mol.line_item;
				if (word.find("gen") != string::npos) ifs >> new_mol.num_generation;
			}

			mol.push_back(new_mol);
		}
	}
	ifs.close();

	// тест 
	for (int i = 0; i < mol.size(); i++) {
		cout << "number of molecule: " << i + 1 << endl;
		cout << mol[i].line_item << endl;
		cout << mol[i].num_generation << endl;
	}
}

void System::ReadMethods() {
	// читаем входной файл - ищем ключевое слово - узнаем тип Геометрии
	ifstream ifs(source, ifstream::in);
	string word;
	while (ifs.good()) {
		ifs >> word;
		if (word.find("method") != string::npos) {
			double tolerance_;
			string name_method;
			// пока не встретим новый блок или не закончится файл
			while (word != "[" && ifs.good()) {
				ifs >> word;
				if (word.find("type") != string::npos) ifs >> name_method;
				if (word.find("tolerance") != string::npos) ifs >> tolerance_;
			}
			//num_iter_, tolerance_, step_
			if (name_method == "gradient") {
				methods.push_back(make_unique<Gradient>(tolerance_));
			}
			else if (name_method == "DFP") {
				methods.push_back(make_unique<DFP>(tolerance_));
			}
			else {
				printf("ERROR: unknown method");
			}

		}
	}
	ifs.close();

	// тест и важный пример того, как использовать 
	// массив с объектами разных методов
	for (BaseOptimTools const &w : methods) {
		cout << w.name << endl;
		cout << w.tolerance << endl;
	}
}





