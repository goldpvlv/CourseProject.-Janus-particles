

//∙ curvature - после этого ключа идет значение кривизны 𝑅
//(используется в torus - е)
//∙ layers_x - (𝑀𝑥)число слоев по 𝑥
//∙ layers_y - (𝑀𝑦)число слоев по 𝑦
//∙ chi_seg - параметр Флори(> 0) для отталкивания между двумя
//молекулами
//∙ num_iter - число шагов в методе сходимости
//∙ tolerance - точность решения
//∙ step - шаг сходимости
//  для молекулы
//∙ ns - длина линейного блока в молекуле
//∙ num_generation - число поколений 𝐺
//∙ sigma - число молекул на элементарную площадь(иногда это называют плотностью прививки 𝜎)
//∙ chi - параметр Флори для взаимодействия с растворителем (х)
//∙ xmin - нижняя граница прививки по 𝑥
//∙ xmax - верхняя граница прививки по 𝑥
//∙ ymin - нижняя граница прививки по 𝑦
//∙ ymax - верхняя граница прививки по 𝑦
//  N - общее число атомов
//  theta - количество вещества = sigma*N



#define _CRT_SECURE_NO_WARNINGS

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


int main() {



	string my_source = "input.txt";
	System system(my_source);

	system.ReadParameters();
	system.SetGeometry();
	system.ReadMolecules();
	system.ReadMethods();
	system.Cycling();


	return 0;
}