#pragma once
#include <fstream>
#include <vector>

namespace element
{
	const int n_func_elem = 9;
	class Element
	{
	public:
		int nodes[4];
		int dof[n_func_elem];
		int ndof;
		int number_of_area;
		int neighbors[4]; //левый, правый, нижний, верхний

		Element& operator=(Element element);
	};

	std::ifstream& operator>>(std::ifstream& is, std::vector <Element>& elements);
}