#pragma once


#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

const double pi = 3.14159265359;
class Geometry {
public:
	mutable int Mx, My;//число слоев по х, y
	double R0, delta_alpha;//кривизна, угол
	void GetValue(int Mx_, int My_, double R0_);//метод для считывания параметров
	Geometry() {};//конструктор класса
	vector < vector < double > >lambda_bb;
	vector < vector < double > >lambda_bn;
	vector < vector < double > >lambda_bf;
	vector < vector < double > >lambda_fb;
	vector < vector < double > >lambda_fn;
	vector < vector < double > >lambda_ff;
	vector < vector < double > >lambda_nb;
	vector < vector < double > >lambda_nn;
	vector < vector < double > >lambda_nf;
	vector < vector < double > >volume;
	vector < vector < double > >square_side;
	vector < vector < double > >square_up;
	vector < vector < double > >square_front;
	vector < vector < double > >square_right;
	vector < vector < double > >square_left;
	virtual void UpdateVolume() const {};      //метод нахождения объема
	virtual void UpdateSquareFront() const {};//метод нахождения площади фронтальной стороны
	virtual void UpdateSquareUp() const {};//метод нахождения площади верхней стороны
	virtual void UpdateSquareLeft() const {};//метод нахождения площади левой стороны
	virtual void UpdateSquareRight() const {};//метод нахождения площади правой стороны
	virtual void UpdateSquareSide() const {};//метод нахождения площади боковой стороны
	virtual void Transposition() const {};//метод нахождения вероятностей перехода
	void AllocateMemory();//инициализация данных
	~Geometry() {};//деконструктор класса
};







class Polar : public Geometry {
public:
	int delta_h = 1;
	void UpdateVolume();
	void UpdateSquareFront();
	void UpdateSquareUp();
	void UpdateSquareSide();
	void Transposition();
};
class Torus : public Geometry {
public:
	void UpdateVolume();
	void UpdateSquareFront();
	void UpdateSquareUp();
	void UpdateSquareLeft();
	void UpdateSquareRight();
	void Transposition();};
