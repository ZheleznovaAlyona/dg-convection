#include "main_solver.h"
#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "element.h"
#include "partition.h"
#include "testing_parameters.h"
#include "parameters.h"
#include "myfunctions.h"
#include "point.h"

using namespace boundary_conditions;
using namespace element;
using namespace partition;
using namespace myvector;
using namespace matrix;
using namespace boundaries;
using namespace logger;
using namespace parameters;
using namespace solver;
using namespace slae;
using namespace testingparameters;
using namespace point;
using namespace basis;

namespace mainsolver
{

	double MainSolver::get_solution_in_point_u(double x, double y, int element_number, MyVector qi)
	{
		//const int nf = 9;
		//int indexes[nf];
		//double u_in_point, qi_local[nf];

		////собираем глобальные номера с элемента
		//for (int j = 0; j < n_func; j++)
		//	indexes[j] = elements[element_number].dof[j];

		////собираем локальный набор весов
		//for (int j = 0; j < n_func; j++)
		//	qi_local[j] = qi[indexes[j]];

		////вычисляем в решение в точке
		//static auto& part = static_cast<Partition>(*this);
		//u_in_point = 0;
		//for (int j = 0; j < n_func; j++)
		//	u_in_point += qi_local[j] * phi_i(j, x, y, element_number, part);

		//return u_in_point;

		//вычисляем в решение в точке
		static auto& part = static_cast<Partition>(*this);
		double u_in_point = 0;
		for (int j = 0; j < elements[element_number].ndof; j++)
			u_in_point += qi[elements[element_number].dof[j]] * phi_i(j, x, y, element_number, part);

		return u_in_point;
	}

	double MainSolver::get_solution_in_point2_u(double x, double y, MyVector qi)
	{
		int element_number = search_element(x, y);
		return get_solution_in_point_u(x, y, element_number, qi);
	}

	void MainSolver::get_vector_solution_in_nodes_u(MyVector qi, MyVector &solution)
	{
		logger.send_message_U();
		//const int nf = 9;
		int indexes[n_func];
		int indexes_nodes[4];
		int size = elements.size();
		double u_local[4], qi_local[n_func];
		double x, y;

		static auto& part = static_cast<Partition>(*this);

		for (int i = 0; i < size; i++)
		{
			//собираем глобальные номера с элемента
			for (int j = 0; j < 4; j++)
				indexes_nodes[j] = elements[i].nodes[j];
			for (int j = 0; j < n_func; j++)
				indexes[j] = elements[i].dof[j];

			//собираем локальный набор весов
			for (int j = 0; j < n_func; j++)
				qi_local[j] = qi[indexes[j]];

			//вычисляем в узлах элемента решение
			for (int j = 0; j < 4; j++)
			{
				u_local[j] = 0;
				x = nodes[indexes_nodes[j]].x;
				y = nodes[indexes_nodes[j]].y;
				for (int k = 0; k < n_func; k++)
					u_local[j] += qi_local[k] * phi_i(k, x, y, i, part);
			}

			//кладём в результирующий вектор
			for (int j = 0; j < 4; j++)
				solution[indexes_nodes[j]] = u_local[j];
		}
	}

	void MainSolver::get_vector_in_nodes_u_an(myvector::MyVector &U_an)
	{
		int indexes_nodes[4];
		int size = elements.size();
		double u_local[4];
		double x, y;

		static auto& part = static_cast<Partition>(*this);

		for (int i = 0; i < size; i++)
		{
			//собираем глобальные номера с элемента
			for (int j = 0; j < 4; j++)
				indexes_nodes[j] = elements[i].nodes[j];

			//вычисляем в узлах элемента решение
			for (int j = 0; j < 4; j++)
			{
				u_local[j] = 0;
				x = nodes[indexes_nodes[j]].x;
				y = nodes[indexes_nodes[j]].y;
				for (int k = 0; k < n_func; k++)
					u_local[j] += calculate_u_analytic(elements[i].number_of_area, x, y);
			}

			//кладём в результирующий вектор
			for (int j = 0; j < 4; j++)
				U_an[indexes_nodes[j]] = u_local[j];
		}
	}
}