#include "partition.h"
#include <iostream>

using namespace std;
using namespace element;
using namespace point;

namespace partition
{
	void Partition::input(ifstream& grid_f_in, ifstream& elements_f_in)
	{
		grid_f_in >> nodes;
		elements_f_in >> elements;
	};

	double Partition::get_hx(int element_number)
	{
		Element element = elements[element_number];
		return nodes[element.nodes[1]].x - nodes[element.nodes[0]].x;
	}

	double Partition::get_hy(int element_number)
	{
		Element element = elements[element_number];
		return nodes[element.nodes[2]].y - nodes[element.nodes[0]].y;
	}

	int Partition::count_unzero_matrix_elements(int slae_size)
	{
		int count = 0;
		for (int i = 0; i < elements.size(); i++)
		{
			int count_local = 0;
			//с собой
			count_local += elements[i].ndof;
			//с соседями
			for (int j = 0; j < 4; j++)
				if (elements[i].neighbors[j] != -1)
					count_local += elements[elements[i].neighbors[j]].ndof;
			count_local *= elements[i].ndof;
			count += count_local;
		}
		count -= slae_size;
		//так как нужно для одного треугольника ввиду симметричности портрета,
		//то необходимо полученное количество поделить на 2
		//std::cout << count / 2;
		return count / 2;
	}

	int Partition::create_unzero_elements_list(int element_number, vector <int> &list)
	{
		int neighbor;

		//свои узлы
		for (int i = 0; i < elements[element_number].ndof; i++)
			list.push_back(elements[element_number].dof[i]);

		//соседей
		for (int j = 0; j < 4; j++)
		{
			neighbor = elements[element_number].neighbors[j];
			if (neighbor != -1)
				for (int i = 0; i < elements[neighbor].ndof; i++)
					list.push_back(elements[neighbor].dof[i]);
		}

		return list.size();
	}

	int Partition::search_element(double x, double y)
	{
		double x_left, x_right, y_low, y_up;
		int size = elements.size();

		for(int i = 0; i < size; i++)
		{
			x_left = nodes[elements[i].nodes[0]].x;
			x_right = nodes[elements[i].nodes[1]].x;
			y_low = nodes[elements[i].nodes[0]].y;
			y_up = nodes[elements[i].nodes[3]].y;
			if(x_left <= x && x <= x_right && y_low <= y && y <= y_up)
				return i;
		}

		return -1;
	}

}