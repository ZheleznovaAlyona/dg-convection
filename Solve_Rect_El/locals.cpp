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

namespace mainsolver
{
	void MainSolver::calculate_locals(int element_number, MyVector q_calc)
	{
		Element element = Partition::elements[element_number];

		calculate_G(element_number);
		calculate_C(element_number, q_calc);
		calculate_M(element_number);
		calculate_F(element_number);
	}

	void MainSolver::calculate_G(int element_number)
	{
		double hx, hy, hx2, hy2, g1, g2;
		Element element = Partition::elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);

		hx = get_hx(element_number);
		hy = get_hy(element_number);
		hx2 = hx * hx;
		hy2 = hy * hy;

		double jacobian = hx * hy / 4.0;

		for (int i = 0; i < element.ndof; i++)
		{
			int id_i = element.dof[i];
			for (int j = 0; j < element.ndof; j++)
			{
				int id_j = element.dof[j];
				g1 = 0;
				g2 = 0;
				for (int k = 0; k < 9; k++)
				{
					double p_ksi = 0.5 + 0.5 * gauss_points[0][k],
						p_etta = 0.5 + 0.5 * gauss_points[1][k];
					g1 += gauss_weights[k] * dphiksi[i](p_ksi, p_etta) *
						dphiksi[j](p_ksi, p_etta);
					g2 += gauss_weights[k] * dphietta[i](p_ksi, p_etta) *
						dphietta[j](p_ksi, p_etta);
				}
				double gij = (g1 * jacobian / hx2 + g2 * jacobian / hy2) * lambda;
				my_slae.A.add_element(id_i, id_j, gij);
				//my_slae.A.mtr[id_i][id_j] += gij;
				//if(i != j) my_slae.A.add_element(id_j, id_i, gij);
			}
		}
	}
	void MainSolver::calculate_C(int element_number, MyVector q_calc)
	{
		double hx, hy;
		Element element = elements[element_number];
		hx = get_hx(element_number);
		hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;
		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		Point a = calculate_a(element.number_of_area);

		for (int i = 0; i < element.ndof; i++)
		{
			int id_i = element.dof[i];
			for (int j = 0; j < element.ndof; j++)
			{
				int id_j = element.dof[j];
				double c1 = 0;
				double c2 = 0;
				for (int k = 0; k < 9; k++)
				{
					double p_ksi = 0.5 + 0.5 * gauss_points[0][k],
						p_etta = 0.5 + 0.5 * gauss_points[1][k];
					c1 += gauss_weights[k] * dphiksi[j](p_ksi, p_etta) *
						phi[i](p_ksi, p_etta);
					c2 += gauss_weights[k] * dphietta[j](p_ksi, p_etta) *
						phi[i](p_ksi, p_etta);
				}
				double cij = c1 * a.x * jacobian / hx + c2 * a.y * jacobian / hy;
				my_slae.A.add_element(id_i, id_j, cij);
				//my_slae.A.mtr[id_i][id_j] += cij;
			}
		}

	}
	void MainSolver::calculate_M(int element_number)
	{
		Element element = elements[element_number];
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;

		double gamma = calculate_gamma(element.number_of_area);

		for (int i = 0; i < element.ndof; i++)
		{
			int id_i = element.dof[i];
			for (int j = 0; j < element.ndof; j++)
			{
				int id_j = element.dof[j];
				double mij = 0;
				for (int k = 0; k < 9; k++)
				{
					double p_ksi = 0.5 + 0.5 * gauss_points[0][k],
						p_etta = 0.5 + 0.5 * gauss_points[1][k];
					mij += gauss_weights[k] * phi[i](p_ksi, p_etta) *
						phi[j](p_ksi, p_etta);
				}
				mij *= jacobian * gamma;
				my_slae.A.add_element(id_i, id_j, mij);
				//my_slae.A.mtr[id_i][id_j] += mij;
			}
		}
	}
	void MainSolver::calculate_F(int element_number)
	{
		double f;
		Element element = elements[element_number];
		int number_of_area = element.number_of_area;
		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;

		for (int i = 0; i < element.ndof; i++)
		{
			int id_i = element.dof[i];
			double fi = 0;
			for (int j = 0; j < 9; j++)
			{
				double p_ksi = 0.5 + 0.5 * gauss_points[0][j],
					p_etta = 0.5 + 0.5 * gauss_points[1][j];
				double p_x = p_ksi * hx + x0, p_y = p_etta * hy + y0;
				f = calculate_f(number_of_area, p_x, p_y);
				fi += f * gauss_weights[j] * phi[i](p_ksi, p_etta);
			}
			fi *= jacobian;
			my_slae.b[id_i] += fi;
		}
	}

}