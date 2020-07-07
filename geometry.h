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
	mutable int Mx, My;
	double R0, alpha, delta_alpha;
	void GetValue(int Mx_, int My_, double R0_);

	Geometry() {};

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

	virtual void UpdateVolume() const {};
	virtual void UpdateSquareFront() const {};
	virtual void UpdateSquareUp() const {};
	virtual void UpdateSquareLeft() const {};
	virtual void UpdateSquareRight() const {};
	virtual void UpdateSquareSide() const {};
	virtual void Transposition() const {};
	void AllocateMemory();
	
	void PrintVolume();
	void Print();
	
	~Geometry() {};
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
	void Transposition();
};