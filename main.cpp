#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include "classes.h"
#include <conio.h>

using namespace std;



int main() {
	string my_source = "input.txt";
	System system(my_source);
	// ������ ������� ���� 
	system.ReadParameters();
	// ������ ���������
	system.SetGeometry();
	// ������ ��������� �������
	system.ReadMolecules();
	// ������ ������ �������
	system.ReadMethods();

	system.GetMethods();
	return 0;
}


