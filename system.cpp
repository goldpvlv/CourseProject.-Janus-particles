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
			}

			new_mol.MyNewMethod();
			mol.push_back(new_mol);
		}
	}
	ifs.close();

	
	/*for (int i = 0; i < mol.size(); i++) {
		cout << "number of molecule: " << i + 1 << endl;
		cout << mol[i].line_item << endl;
		cout << mol[i].num_generation << endl;
	}*/
};

void System::ReadMethods() {

	ifstream ifs(source, ifstream::in);
	string word;
	while (ifs.good()) {
		ifs >> word;
		if (word.find("method") != string::npos) {
			double tolerance_, nu_, _num_iter;
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
				methods.push_back(make_unique<Gradient>(nu_));
				methods.push_back(make_unique<Gradient>(_num_iter));
			}
			else if (name_method == "DFP") {
				methods.push_back(make_unique<DFP>(tolerance_));
				methods.push_back(make_unique<DFP>(nu_));
				methods.push_back(make_unique<DFP>(_num_iter));
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
		cout << w.nu << endl;
		cout << w.num_iter << endl;
	}
}



void System::Cycling() {
	for (BaseOptimTools const &scheme : methods) {
		cout << scheme.name << endl;

		double length_of_grad = 0.0;
		int step = 0.0;

		Function();
		length_of_grad = scheme.SetGradFirst(grad);

		do {
			scheme.UpdateX(u, grad);
			Function();
			length_of_grad = scheme.SetGradRegular(grad);

			step++;
			cout << step << "     " << length_of_grad << endl;
		} while ((length_of_grad > scheme.tolerance) && (step < scheme.num_iter));
	}
};