#pragma once


#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <stdio.h>


#include "geometry.h"


using namespace std;
using vector3d = vector<vector<vector<double>>>;

class Molecule {
	friend class System;
public:
	double theta, sigma, chi;
	int layers_x, layers_y, xmin, xmax, ymin, ymax, num_generation, num_atoms, ns;
	
	void SetParameters();//метод считывания параметров молекулы
	void AllocateMemory(int layers_x_, int layers_y_, int M);//метод инициализации данных
	double q;
	vector3d Gback;
	vector3d Gforw;
	vector <double> u;
	vector <vector <double>> G;
	vector <double> fi_side;
	vector<vector<double>> fi_p;
	vector<vector<double>> multipliers;

	void FindG();//метод для нахождения распределения Больцмана
	void FindGforw(Geometry geo);//метод нахождения прямого пропагатора
	void FindGback(Geometry geo);//метод нахождения обратного пропагатора
	void FindQ(Geometry geo);//метод нахождения статической суммы
	void FindFiP();//метод нахождения профилей плотности молекулы
	void FindFiSide(Geometry geo);//метод нахождения профиля плотности окружения

};
