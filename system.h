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
const double T = 298.15;
const double k_B = 1.38065812e-23;
const double eps = 80;
const double eps_0 = 8.85418785e-12;
const double e = 1.6021773349e-19;
const double b = 3.1e-10;
const double pi = 3.14159265359;

template<typename T> struct Box : std::unique_ptr<T> {
	using std::unique_ptr<T>::unique_ptr;
	operator T &() const { return **this; }};
class System
{ 
	
public:
	
	double f_n = e * e / pi*eps * eps_0*b*k_B*T;
	string source;
	vector <double> u;
	vector <double> grad;

	void ReadParameters();//метод для считывания параметров системы

	int layers_x;
	int layers_y;
	int M;
	double Rcurv, chi_seg, branch;
	string geometry_name;
	System(string _source);//конструктор класса
	Geometry geo;//объявление объекта типа геометрии
	void SetGeometry();//метод для считывания параметров геометрии
	vector<Molecule> mol;//массив молекул
	void ReadMolecules();//метод для считывания параметров молекул
	vector<Box<BaseOptimTools>> methods; //шаблон для реализации методов
	void ReadMethods();//метод для считывания 
	int get_branch() const;
	void Function();// метод самосогласованного поля
	void Cycling();// метод поиска градиента
	void Output();//метод вывода в файл

private:

	void AllocateMemory();//метод инициализация данных
	void set_fi_salt(const Molecule &fi_salt);
	void FindVolumeRatioSolv();//подсчет объемной доли растворителя
	void FindFiSolv();//метод расчета профиля плотности растворителя
	void FindFiTotal();//метод расчета общего профиля плотности
	void FindFiSide();//метод расчета профиля плотности окружности для растворителя
	void FindLagrangeMultipliers(int t);//метод нахождения множителей Лагранжа для молекул
	void FindLagrangeMultipliers();//метод нахождения множителей Лагранжа для растворителя
	void FindGrad();//метод расчета градиента
	void FindCharge(Geometry geo);// расчет заряда
	void ElectrostaticPotential();//значение электростатического потенциала
	void PoissonBoltzmannEquation();//уравнение Пуассона-Больцмана

	double _fi_salt;
	vector <double> fi_side;
	vector <vector <double>> ksi_0;
	vector <vector <double>> ksi;
	vector <vector <double>> q;
	vector < vector < double > > fi_solv;
	vector < vector < double > > fi_total;
	vector <vector<double>> multipliers; 
	vector <vector<double>> middle_multipliers;

};
