#include "boundaries.h"
#include "element.h"
#include "myfunctions.h"

using namespace std;
using namespace element;
using namespace partition;
using namespace parameters;
using namespace basis;
using namespace matrix;

namespace boundaries
{

	void OuterBoundaries::initialize(Formulation formulation)
	{
		mu = 0;
		form = formulation;
	}

	void OuterBoundaries::calculate_outer_boundaries(int element_number, matrix::Matrix& A)
	{
		Element element = elements[element_number];

		int last_node = nodes.size() - 1;
		int left_low_corner_node = element.nodes[0];
		int right_up_corner_node = element.nodes[3];

		if (nodes[left_low_corner_node].x == nodes[0].x)
		{
			calculate_ES_out_left(element_number, A);
		}//вертикальная левая граница
		if (nodes[right_up_corner_node].x == nodes[last_node].x)
		{
			calculate_ES_out_right(element_number, A);
		}//вертикальная правая граница
		if (nodes[left_low_corner_node].y == nodes[0].y)
		{
			calculate_ES_out_low(element_number, A);
		}//горизонтальная нижняя граница
		if (nodes[right_up_corner_node].y == nodes[last_node].y)
		{
			calculate_ES_out_up(element_number, A);
		}//горизонтальная верхняя граница
	}

	void OuterBoundaries::calculate_ES_out_left(int element_number, Matrix& A)
	{
		//const int nf = 9;
		double AK[n_func][n_func], S[n_func][n_func];
		for (int i = 0; i < n_func; i++)
		{
			memset(&AK[i][0], 0, sizeof(double) * n_func);
			memset(&S[i][0], 0, sizeof(double) * n_func);
		}

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);
		double a = hy / hx;
		a *= lambda * 0.5;//якобиан*lambda *dksi/dx
		int id_i, id_j;

		mu = 0;
		double jacobian = 0.5 * hy;
		double st = jacobian * mu;

		if (form == BaumannOden || form == NIPG)
		{
			for (int i = 0; i < n_func; i++)
			{
				for (int j = 0; j < n_func; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						double p_etta = 0.5 + 0.5 * gauss_points_1[k];
						AK[i][j] += gauss_weights_1[k] * (phi[j](0, p_etta) *
							dphiksi[i](0, p_etta) - phi[i](0, p_etta) *
							dphiksi[j](0, p_etta));
						S[i][j] += gauss_weights_1[k] * phi[j](0, p_etta) * phi[i](0, p_etta);
					}
					AK[i][j] *= -a;
					S[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
			{
				id_i = element.dof[i];
				for (int j = 0; j < n_func; j++)
				{
					id_j = element.dof[j];
					A.add_element(id_i, id_j, AK[i][j] + S[i][j]);
					//A.mtr[id_i][id_j] += AK[i][j] + S[i][j];
				}
			}
		}

		if (form == IP)
		{
			for (int i = 0; i < n_func; i++)
			{
				for (int j = 0; j < n_func; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						double p_etta = 0.5 + 0.5 * gauss_points_1[k];
						AK[i][j] += gauss_weights_1[k] * (phi[j](0, p_etta) *
							dphiksi[i](0, p_etta) + phi[i](0, p_etta) *
							dphiksi[j](0, p_etta));
						S[i][j] += gauss_weights_1[k] * phi[j](0, p_etta) * phi[i](0, p_etta);
					}
					AK[i][j] *= -a;
					S[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
			{
				id_i = element.dof[i];
				for (int j = 0; j < n_func; j++)
				{
					id_j = element.dof[j];
					A.add_element(id_i, id_j, -AK[i][j] + S[i][j]);
					//A.mtr[id_i][id_j] += -AK[i][j] + S[i][j];
				}
			}
		}

	}
	void OuterBoundaries::calculate_ES_out_right(int element_number, Matrix& A)
	{
		//const int nf = 9;
		double AN[n_func][n_func], S[n_func][n_func];
		for (int i = 0; i < n_func; i++)
		{
			memset(&AN[i][0], 0, sizeof(double) * n_func);
			memset(&S[i][0], 0, sizeof(double) * n_func);
		}

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);
		double a = hy / hx;
		a *= lambda * 0.5;//якобиан*lambda *dksi/dx
		int id_i, id_j;

		mu = 0;
		double jacobian = 0.5 * hy;
		double st = jacobian * mu;

		if (form == BaumannOden || form == NIPG)
		{
			for (int i = 0; i < n_func; i++)
			{
				for (int j = 0; j < n_func; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						double p_etta = 0.5 + 0.5 * gauss_points_1[k];
						AN[i][j] += gauss_weights_1[k] * (phi[j](1, p_etta) *
							dphiksi[i](1, p_etta) - phi[i](1, p_etta) *
							dphiksi[j](1, p_etta));
						S[i][j] += gauss_weights_1[k] * phi[j](1, p_etta) * phi[i](1, p_etta);
					}
					AN[i][j] *= a;
					S[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
			{
				id_i = element.dof[i];
				for (int j = 0; j < n_func; j++)
				{
					id_j = element.dof[j];
					A.add_element(id_i, id_j, AN[i][j] + S[i][j]);
					//A.mtr[id_i][id_j] += AN[i][j] + S[i][j];
				}
			}
		}

		if (form == IP)
		{
			for (int i = 0; i < n_func; i++)
			{
				for (int j = 0; j < n_func; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						double p_etta = 0.5 + 0.5 * gauss_points_1[k];
						AN[i][j] += gauss_weights_1[k] * (phi[j](1, p_etta) *
							dphiksi[i](1, p_etta) + phi[i](1, p_etta) *
							dphiksi[j](1, p_etta));
						S[i][j] += gauss_weights_1[k] * phi[j](1, p_etta) * phi[i](1, p_etta);
					}
					AN[i][j] *= a;
					S[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
			{
				id_i = element.dof[i];
				for (int j = 0; j < n_func; j++)
				{
					id_j = element.dof[j];
					A.add_element(id_i, id_j, -AN[i][j] + S[i][j]);
					//A.mtr[id_i][id_j] += -AN[i][j] + S[i][j];
				}
			}
		}
	}
	void OuterBoundaries::calculate_ES_out_low(int element_number, Matrix& A)
	{
		//const int nf = 9;
		double AK[n_func][n_func], S[n_func][n_func];
		for (int i = 0; i < n_func; i++)
		{
			memset(&AK[i][0], 0, sizeof(double) * n_func);
			memset(&S[i][0], 0, sizeof(double) * n_func);
		}

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);
		double a = hx / hy;
		a *= lambda * 0.5;//якобиан*lambda *detta/dy
		int id_i, id_j;

		mu = 0;
		double jacobian = 0.5 * hx;
		double st = jacobian * mu;

		if (form == BaumannOden || form == NIPG)
		{
			for (int i = 0; i < n_func; i++)
			{
				for (int j = 0; j < n_func; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						double p_ksi = 0.5 + 0.5 * gauss_points_1[k];
						AK[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 0) *
							dphietta[i](p_ksi, 0) - phi[i](p_ksi, 0) *
							dphietta[j](p_ksi, 0));
						S[i][j] += gauss_weights_1[k] * phi[j](p_ksi, 0) * phi[i](p_ksi, 0);
					}
					AK[i][j] *= -a;
					S[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
			{
				id_i = element.dof[i];
				for (int j = 0; j < n_func; j++)
				{
					id_j = element.dof[j];
					A.add_element(id_i, id_j, AK[i][j] + S[i][j]);
					//A.mtr[id_i][id_j] += AK[i][j] + S[i][j];
				}
			}
		}

		if (form == IP)
		{
			for (int i = 0; i < n_func; i++)
			{
				for (int j = 0; j < n_func; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						double p_ksi = 0.5 + 0.5 * gauss_points_1[k];
						AK[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 0) *
							dphietta[i](p_ksi, 0) + phi[i](p_ksi, 0) *
							dphietta[j](p_ksi, 0));
						S[i][j] += gauss_weights_1[k] * phi[j](p_ksi, 0) * phi[i](p_ksi, 0);
					}
					AK[i][j] *= -a;
					S[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
			{
				id_i = element.dof[i];
				for (int j = 0; j < n_func; j++)
				{
					id_j = element.dof[j];
					A.add_element(id_i, id_j, -AK[i][j] + S[i][j]);
					//A.mtr[id_i][id_j] += -AK[i][j] + S[i][j];
				}
			}
		}
	}
	void OuterBoundaries::calculate_ES_out_up(int element_number, Matrix& A)
	{
		//const int nf = 9;
		double AN[n_func][n_func], S[n_func][n_func];
		for (int i = 0; i < n_func; i++)
		{
			memset(&AN[i][0], 0, sizeof(double) * n_func);
			memset(&S[i][0], 0, sizeof(double) * n_func);
		}

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);
		double a = hx / hy;
		a *= lambda * 0.5;//якобиан*lambda *detta/dy
		int id_i, id_j;

		mu = 0;
		double jacobian = 0.5 * hx;
		double st = jacobian * mu;

		if (form == BaumannOden || form == NIPG)
		{
			for (int i = 0; i < n_func; i++)
			{
				for (int j = 0; j < n_func; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						double p_ksi = 0.5 + 0.5 * gauss_points_1[k];
						AN[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 1) *
							dphietta[i](p_ksi, 1) - phi[i](p_ksi, 1) *
							dphietta[j](p_ksi, 1));
						S[i][j] += gauss_weights_1[k] * phi[j](p_ksi, 1) * phi[i](p_ksi, 1);
					}
					AN[i][j] *= a;
					S[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
			{
				id_i = element.dof[i];
				for (int j = 0; j < n_func; j++)
				{
					id_j = element.dof[j];
					A.add_element(id_i, id_j, AN[i][j] + S[i][j]);
					//A.mtr[id_i][id_j] += AN[i][j] + S[i][j];
				}
			}
		}

		if (form == IP)
		{
			for (int i = 0; i < n_func; i++)
			{
				for (int j = 0; j < n_func; j++)
				{
					for (int k = 0; k < 3; k++)
					{
						double p_ksi = 0.5 + 0.5 * gauss_points_1[k];
						AN[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 1) *
							dphietta[i](p_ksi, 1) + phi[i](p_ksi, 1) *
							dphietta[j](p_ksi, 1));
						S[i][j] += gauss_weights_1[k] * phi[j](p_ksi, 1) * phi[i](p_ksi, 1);
					}
					AN[i][j] *= a;
					S[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
			{
				id_i = element.dof[i];
				for (int j = 0; j < n_func; j++)
				{
					id_j = element.dof[j];
					A.add_element(id_i, id_j, -AN[i][j] + S[i][j]);
					//A.mtr[id_i][id_j] += -AN[i][j] + S[i][j];
				}
			}
		}
	}
}