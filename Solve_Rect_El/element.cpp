#include "element.h"
#include "basis.h"

using namespace std;
using namespace basis;

namespace element
{
	Element& Element::operator=(Element element)
	{
		ndof = element.ndof;
		for(int i = 0; i < 4; i++)
			nodes[i] = element.nodes[i];
		for(int i = 0; i < ndof; i++)
			dof[i] = element.dof[i];
		number_of_area = element.number_of_area;
		for(int i = 0; i < 4; i++)
		neighbors[i] = element.neighbors[i];

		return *this;
	}

	ifstream& operator>>(ifstream& is, vector <Element>& elements)
	{
		int tmp;
		Element element_tmp;

		is >> tmp;
		elements.reserve(tmp);

		for(int i = 0; i < tmp; i++)
		{
			is >> element_tmp.number_of_area;
			for(int j = 0; j < 4; j++)
				is >> element_tmp.nodes[j];
			for(int j = 0; j < 4; j++)
				is >> element_tmp.neighbors[j];

			element_tmp.ndof = n_func;
			element_tmp.dof[0] = element_tmp.ndof * i;
			for (int j = 1; j < element_tmp.ndof; j++)
				element_tmp.dof[j] = element_tmp.dof[0] + j;
			elements.push_back(element_tmp);
		}

		//int count = -1;
		//for (int j = 0; j < 16; j++)
		//{
		//	for (int i = 0; i < 32; i++)
		//	{

		//		count++;
		//		for (int iy = 0; iy < 3; iy++)
		//		{
		//			int jj = 3 * j + iy;
		//			for (int ix = 0; ix < 3; ix++)
		//			{
		//				int ii = 3 * i + ix;
		//				elements[count].dof[iy * 3 + ix] = 3 * jj * 32 + ii;
		//			}
		//		}
		//	}
		//}

		return is;
	}
}