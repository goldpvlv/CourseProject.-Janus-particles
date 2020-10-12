#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <stdio.h>

#include "geometry.h"
#include "molecule.h"
#include "method.h"   

using namespace std;
using vector3d = vector<vector<vector<double>>>;

template<typename T> struct Box : std::unique_ptr<T> {
	using std::unique_ptr<T>::unique_ptr;
	operator T &() const { return **this; }};
class System
{ 
public:
	
	string source;
	vector <double> u;
	vector <double> grad;

	void ReadParameters();//метод для считывания параметров системы

	int layers_x;
	int layers_y;
	int M;
	double Rcurv, chi_seg;
	string geometry_name;
	System(string _source);//конструктор класса
	Geometry geo;//объявление объекта типа геометрии
	void SetGeometry();//метод для считывания параметров геометрии
	vector<Molecule> mol;//массив молекул
	void ReadMolecules();//метод для считывания параметров молекул
	vector<Box<BaseOptimTools>> methods; //шаблон для реализации методов
	void ReadMethods();//метод для считывания 
	void Function();// метод самосогласованного поля
	void Cycling();// метод поиска градиента
	void Output();//метод вывода в файл

private:

	void AllocateMemory();//метод инициализация данных

	void FindFiSolv();//метод расчета профиля плотности растворителя
	void FindFiTotal();//метод расчета общего профиля плотности
	void FindFiSide();//метод расчета профиля плотности окружности для растворителя
	void FindLagrangeMultipliers(int t);//метод нахождения множителей Лагранжа для молекул
	void FindLagrangeMultipliers();//метод нахождения множителей Лагранжа для растворителя
	void FindGrad();//метод расчета градиента

	vector <double> fi_side;
	vector < vector < double > > fi_solv;
	vector < vector < double > > fi_total;
	vector <vector<double>> multipliers; 
	vector <vector<double>> middle_multipliers;

};
