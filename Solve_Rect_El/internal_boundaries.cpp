#include "boundaries.h"
#include "element.h"
#include "myfunctions.h"
#include <iostream>

using namespace std;
using namespace element;
using namespace partition;
using namespace parameters;
using namespace basis;
using namespace matrix;

namespace boundaries
{
	void InternalBoundaries::initialize(Formulation formulation)
	{
		mu = 0;
		form = formulation;
	}

	void InternalBoundaries::calculate_internal_boundaries(int element_number, Matrix& A)
	{
		Element element = elements[element_number];
		int neighbor_element;

		for(int k = 0; k < 4; k++)
		{
			neighbor_element = element.neighbors[k];
			//существующий соседний элемент с бОльшим номером
			if(neighbor_element > element_number)
			{
				//левый/правый сосед->вертикальная граница
				if(k == 0 || k == 1)
				{
					calculate_E_vertical(element_number, neighbor_element, A);
				}
				else
				{
					calculate_E_horizontal(element_number, neighbor_element, A);
				}
			}
		}
	}

	void InternalBoundaries::calculate_E_horizontal(int element_number1, int element_number2, matrix::Matrix& A)
	{
		//const int nf = 9;
		double AK[n_func][n_func], AN[n_func][n_func], BK[n_func][n_func], BN[n_func][n_func],
			   SNN[n_func][n_func], SNK[n_func][n_func], SKN[n_func][n_func], SKK[n_func][n_func];
		vector <vector<double>> ES;
		ES.resize(n_func * 2);

		for (int i = 0; i < n_func; i++)
		{
			memset(&AN[i][0], 0, sizeof(double) * n_func);
			memset(&AK[i][0], 0, sizeof(double) * n_func);
			memset(&BN[i][0], 0, sizeof(double) * n_func);
			memset(&BK[i][0], 0, sizeof(double) * n_func);
			memset(&SNN[i][0], 0, sizeof(double) * n_func);
			memset(&SNK[i][0], 0, sizeof(double) * n_func);
			memset(&SKK[i][0], 0, sizeof(double) * n_func);
			memset(&SKN[i][0], 0, sizeof(double) * n_func);
			initialize_vector(ES[i], n_func * 2);
			initialize_vector(ES[i + n_func], n_func * 2);
		}
		

		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double lambda = calculate_lambda(element.number_of_area);
		double lambda_2 = calculate_lambda(element_2.number_of_area);
		double hx = get_hx(element_number1);
		double hy = get_hy(element_number1);
		double hy_2 = get_hy(element_number2);
		double a = hx / hy;
		double a_lambda = a * lambda * 0.25;//якобиан*0.5*lambda *detta/dy
		double b = hx / hy_2;
		double b_lambda = b * lambda_2 * 0.25;
		double aa_lambda = 0.25 * hx * lambda; //якобиан*0.5*lambda
		double bb_lambda = 0.25 * hx * lambda_2;

		if(form == NIPG || form == IP) 
		{ 
			double na = calculate_a(element.number_of_area).norm();
			if(fabs(na) > 1e-10) mu = 4 * 2 * lambda  *  na / hx;
			else mu = 4 * 2 * lambda / hx;
		}	
		else 
			if (form == BaumannOden)
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
						AK[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 0) *
							dphietta[i](p_ksi, 0) - phi[i](p_ksi, 0) *
							dphietta[j](p_ksi, 0));
						BN[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 0) *
							dphietta[i](p_ksi, 1) / hy + phi[i](p_ksi, 1) *
							dphietta[j](p_ksi, 0) / hy_2);
						BK[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 1) *
							dphietta[i](p_ksi, 0) / hy_2 + phi[i](p_ksi, 0) *
							dphietta[j](p_ksi, 1) / hy);
						SNN[i][j] += gauss_weights_1[k] * phi[j](p_ksi, 1) * phi[i](p_ksi, 1);
						SNK[i][j] -= gauss_weights_1[k] * phi[j](p_ksi, 0) * phi[i](p_ksi, 1);
						SKK[i][j] += gauss_weights_1[k] * phi[j](p_ksi, 0) * phi[i](p_ksi, 0);
						SKN[i][j] -= gauss_weights_1[k] * phi[j](p_ksi, 1) * phi[i](p_ksi, 0);

					}
					AN[i][j] *= a_lambda;
					AK[i][j] *= -b_lambda;
					BN[i][j] *= -aa_lambda;
					BK[i][j] *= bb_lambda;
					SNN[i][j] *= st;
					SNK[i][j] *= st;
					SKK[i][j] *= st;
					SKN[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
				for (int j = 0; j < n_func; j++)
				{
					ES[i][j] = AN[i][j] + SNN[i][j];
					ES[i + n_func][j + n_func] = AK[i][j] + SKK[i][j];
					ES[i][j + n_func] = BN[i][j] + SNK[i][j];
					ES[i + n_func][j] = BK[i][j] + SKN[i][j];
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
						AK[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 0) *
							dphietta[i](p_ksi, 0) + phi[i](p_ksi, 0) *
							dphietta[j](p_ksi, 0));
						BN[i][j] += gauss_weights_1[k] * (-phi[j](p_ksi, 0) *
							dphietta[i](p_ksi, 1) / hy + phi[i](p_ksi, 1) *
							dphietta[j](p_ksi, 0) / hy_2);
						BK[i][j] += gauss_weights_1[k] * (phi[j](p_ksi, 1) *
							dphietta[i](p_ksi, 0) / hy_2 - phi[i](p_ksi, 0) *
							dphietta[j](p_ksi, 1) / hy);
						SNN[i][j] += gauss_weights_1[k] * phi[j](p_ksi, 1) * phi[i](p_ksi, 1);
						SNK[i][j] -= gauss_weights_1[k] * phi[j](p_ksi, 0) * phi[i](p_ksi, 1);
						SKK[i][j] += gauss_weights_1[k] * phi[j](p_ksi, 0) * phi[i](p_ksi, 0);
						SKN[i][j] -= gauss_weights_1[k] * phi[j](p_ksi, 1) * phi[i](p_ksi, 0);

					}
					AN[i][j] *= a_lambda;
					AK[i][j] *= -b_lambda;
					BN[i][j] *= aa_lambda;
					BK[i][j] *= bb_lambda;
					SNN[i][j] *= st;
					SNK[i][j] *= st;
					SKK[i][j] *= st;
					SKN[i][j] *= st;
				}
			}

			for (int i = 0; i < n_func; i++)
				for (int j = 0; j < n_func; j++)
				{
					ES[i][j] = -AN[i][j] + SNN[i][j];
					ES[i + n_func][j + n_func] = -AK[i][j] + SKK[i][j];
					ES[i][j + n_func] = -BN[i][j] + SNK[i][j];
					ES[i + n_func][j] = -BK[i][j] + SKN[i][j];
				}
		}

		add_E_to_global(element_number1, element_number2, A, ES);
	}

	void InternalBoundaries::calculate_E_vertical(int element_number1, int element_number2, matrix::Matrix& A)
	{
		//const int nf = 9;
		double AK[n_func][n_func], AN[n_func][n_func], BK[n_func][n_func], BN[n_func][n_func],
			SNN[n_func][n_func], SNK[n_func][n_func], SKN[n_func][n_func], SKK[n_func][n_func];
		vector <vector<double>> ES;
		ES.resize(n_func * 2);

		for (int i = 0; i < n_func; i++)
		{
			memset(&AN[i][0], 0, sizeof(double) * n_func);
			memset(&AK[i][0], 0, sizeof(double) * n_func);
			memset(&BN[i][0], 0, sizeof(double) * n_func);
			memset(&BK[i][0], 0, sizeof(double) * n_func);
			memset(&SNN[i][0], 0, sizeof(double) * n_func);
			memset(&SNK[i][0], 0, sizeof(double) * n_func);
			memset(&SKK[i][0], 0, sizeof(double) * n_func);
			memset(&SKN[i][0], 0, sizeof(double) * n_func);
			initialize_vector(ES[i], n_func * 2);
			initialize_vector(ES[i + n_func], n_func * 2);
		}

		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double lambda = calculate_lambda(element.number_of_area);
		double lambda_2 = calculate_lambda(element_2.number_of_area);
		double hx = get_hx(element_number1);
		double hy = get_hy(element_number1);
		double hx_2 = get_hx(element_number2);
		double a = hy / hx;
		double a_lambda = a * lambda * 0.25; //якобиан*0.5*lambda *dksi/dx
		double b = hy / hx_2;
		double b_lambda = b * lambda_2 * 0.25;
		double aa_lambda = 0.25 * hy * lambda; //якобиан*0.5*lambda
		double bb_lambda = 0.25 * hy * lambda_2;

		if (form == NIPG || form == IP)
		{
			double na = calculate_a(element.number_of_area).norm();
			if (fabs(na) > 1e-10) mu = 4 * 2 * lambda  *  na / hy;
			else mu = 4 * 2 * lambda / hy;
		}
		else
			if (form == BaumannOden)
				mu = 0;
		double jacobian = 0.5 * hy;
		double st = jacobian * mu;

		if(form == BaumannOden || form == NIPG)
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
						AK[i][j] += gauss_weights_1[k] * (phi[j](0, p_etta) *
									dphiksi[i](0, p_etta) - phi[i](0, p_etta) *
									dphiksi[j](0, p_etta));
						BN[i][j] += gauss_weights_1[k] * (phi[j](0, p_etta) *
									dphiksi[i](1, p_etta) / hx + phi[i](1, p_etta) *
									dphiksi[j](0, p_etta) / hx_2);
						BK[i][j] += gauss_weights_1[k] * (phi[j](1, p_etta) *
									dphiksi[i](0, p_etta) / hx_2 + phi[i](0, p_etta) *
									dphiksi[j](1, p_etta) / hx);
						SNN[i][j] += gauss_weights_1[k] * phi[j](1, p_etta) * phi[i](1, p_etta);
						SNK[i][j] -= gauss_weights_1[k] * phi[j](0, p_etta) * phi[i](1, p_etta);
						SKK[i][j] += gauss_weights_1[k] * phi[j](0, p_etta) * phi[i](0, p_etta);
						SKN[i][j] -= gauss_weights_1[k] * phi[j](1, p_etta) * phi[i](0, p_etta);
					}
					AN[i][j] *= a_lambda;
					AK[i][j] *= -b_lambda;
					BN[i][j] *= -aa_lambda;
					BK[i][j] *= bb_lambda;
					SNN[i][j] *= st;
					SNK[i][j] *= st;
					SKK[i][j] *= st;
					SKN[i][j] *= st;
				}
			}
			for (int i = 0; i < n_func; i++)
				for (int j = 0; j < n_func; j++)
				{
					ES[i][j] = AN[i][j] + SNN[i][j];
					ES[i + n_func][j + n_func] = AK[i][j] + SKK[i][j];
					ES[i][j + n_func] = BN[i][j] + SNK[i][j];
					ES[i + n_func][j] = BK[i][j] + SKN[i][j];
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
						AK[i][j] += gauss_weights_1[k] * (phi[j](0, p_etta) *
							dphiksi[i](0, p_etta) + phi[i](0, p_etta) *
							dphiksi[j](0, p_etta));
						BN[i][j] += gauss_weights_1[k] * (-phi[j](0, p_etta) *
							dphiksi[i](1, p_etta) / hx + phi[i](1, p_etta) *
							dphiksi[j](0, p_etta) / hx_2);
						BK[i][j] += gauss_weights_1[k] * (phi[j](1, p_etta) *
							dphiksi[i](0, p_etta) / hx_2 - phi[i](0, p_etta) *
							dphiksi[j](1, p_etta) / hx);
						SNN[i][j] += gauss_weights_1[k] * phi[j](1, p_etta) * phi[i](1, p_etta);
						SNK[i][j] -= gauss_weights_1[k] * phi[j](0, p_etta) * phi[i](1, p_etta);
						SKK[i][j] += gauss_weights_1[k] * phi[j](0, p_etta) * phi[i](0, p_etta);
						SKN[i][j] -= gauss_weights_1[k] * phi[j](1, p_etta) * phi[i](0, p_etta);
					}
					AN[i][j] *= a_lambda;
					AK[i][j] *= -b_lambda;
					BN[i][j] *= aa_lambda;
					BK[i][j] *= bb_lambda;
					SNN[i][j] *= st;
					SNK[i][j] *= st;
					SKK[i][j] *= st;
					SKN[i][j] *= st;
				}
			}
			for (int i = 0; i < n_func; i++)
				for (int j = 0; j < n_func; j++)
				{
					ES[i][j] = -AN[i][j] + SNN[i][j];
					ES[i + n_func][j + n_func] = -AK[i][j] + SKK[i][j];
					ES[i][j + n_func] = -BN[i][j] + SNK[i][j];
					ES[i + n_func][j] = -BK[i][j] + SKN[i][j];
				}
		}

		add_E_to_global(element_number1, element_number2, A, ES);
	}
	void InternalBoundaries::add_E_to_global(int element_number, int neighbor_element_number, Matrix& A, vector<vector<double>>& E)
	{
		int id_i, id_j;
		Element element = elements[element_number];
		Element neighbor_element = elements[neighbor_element_number];

		for(int i = 0; i < n_func; i++)
		{
			id_i = element.dof[i];
			for(int j = 0; j < n_func; j++)
			{				
				id_j = element.dof[j];
				A.add_element(id_i, id_j, E[i][j]); 
				//A.mtr[id_i][id_j] += E[i][j];
			}

			for(int j = n_func; j < n_func * 2; j++)
			{				
				id_j = neighbor_element.dof[j - n_func];
				A.add_element(id_i, id_j, E[i][j]); 
				//A.mtr[id_i][id_j] += E[i][j];
			}
		}

		for(int i = n_func; i < n_func * 2; i++)
		{
			id_i = neighbor_element.dof[i - n_func];
			for(int j = 0; j < n_func; j++)
			{				
				id_j = element.dof[j];
				A.add_element(id_i, id_j, E[i][j]); 
				//A.mtr[id_i][id_j] += E[i][j];
			}

			for(int j = n_func; j < n_func * 2; j++)
			{				
				id_j = neighbor_element.dof[j - n_func];
				A.add_element(id_i, id_j, E[i][j]); 
				//A.mtr[id_i][id_j] += E[i][j];
			}
		}
	}
}


