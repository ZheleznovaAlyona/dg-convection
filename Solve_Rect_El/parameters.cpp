#include "parameters.h"
#include "testing_parameters.h"
#include "myfunctions.h"

using namespace testingparameters;
using namespace point;

namespace parameters
{
	double Parameters::Ug(int formula_number, double x, double y)
	{
		if (Testing_parameters::test == 1)
			switch (formula_number)
			{
			case 0: return x + y; break;
			case 1:	return x + y; break;
			}

		if (Testing_parameters::test == 2)
			switch (formula_number)
			{
			case 0: return y * y + x * x; break;
			case 1:	return y * y + x * x; break;
			case 2:	return y * y + x * x; break;
			case 3:	return y * y + x * x; break;
			}

		if (Testing_parameters::test == 3)
			switch (formula_number)
			{
			case 0: return 0.0; break;
			case 1:	return 0.0; break;
			case 2:	return 0.0; break;
			case 3:	return 0.0; break;
			}

		if (Testing_parameters::test == 4)
			switch (formula_number)
			{
			case 0: return x * x * y * y; break;
			case 1:	return x * x * y * y; break;
			}

		if (Testing_parameters::test == 5)
			switch (formula_number)
			{
			case 0: return 0.0; break;
			case 1:	return 0.0; break;
			}

		if (Testing_parameters::test == 6)
			switch (formula_number)
			{
			case 0: return 0.0; break;
			case 1:	return 1.0; break;
			}

		if (Testing_parameters::test == 7)
			switch (formula_number)
			{
			case 0: return 0.0; break;
			case 1:	return 1.0; break;
			case 2:	return 2 * y; break;
			}

		if (Testing_parameters::test == 8)
			switch (formula_number)
			{
			case 0: return 0.0; break;
			case 1:	return 1.0; break;
			}

		return 1.0;
	}

	double Parameters::calculate_f(int area_number, double x, double y)
	{
		Point a = calculate_a(area_number);
		if (Testing_parameters::test == 1)
			switch (area_number)
			{
			case 0: return 0.0; break;
			case 1:	return 0.0; break;
			}

		if (Testing_parameters::test == 2)
			switch (area_number)
			{
			case 0: return -4.0; break;
			case 1:	return -4.0; break;
			}

		if (Testing_parameters::test == 3)
			switch (area_number)
			{
			case 0: return 2 * (x + y) - 2 * (x * x + y * y); break;
			case 1:	return 2 * (x + y) - 2 * (x * x + y * y); break;
			}

		if (Testing_parameters::test == 4)
			switch (area_number)
			{
			case 0: return -2.0 * calculate_lambda(0) * (x * x + y * y) + a.x * 2 * x * y * y + a.y * 2 * y * x * x; break;
			case 1:	return -2.0 * calculate_lambda(1) * (x * x + y * y) + a.x * 2 * x * y * y + a.y * 2 * y * x * x; break;
			}

		if (Testing_parameters::test == 5)
		{
			double lambda = calculate_lambda(area_number);
			double sqrt_lambda = sqrt(lambda);
			double u_tmp = 0.25 * (sin(4 * M_PI * x) + 2) * (1 - exp(-x / sqrt_lambda)) *
				(1 - exp((x - 1) / sqrt_lambda)) * (1 - exp(-y / sqrt_lambda))
				* (1 - exp((-y - 1) / sqrt_lambda));


			switch (area_number)
			{
			case 0: return -lambda * (-4 * sin(4 * M_PI * x) * M_PI * M_PI *
				(1 - exp(-x / sqrt_lambda)) * (1 - exp((x - 1) / sqrt_lambda)) *
				(1 - exp(-y / sqrt_lambda)) * (1 - exp((-y - 1.0) / sqrt_lambda)) *
				(3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) + 2.0 *
				cos(4.0 * M_PI * x) * M_PI / sqrt_lambda * exp(-x / sqrt_lambda) *
				(1.0 - exp((x - 1.0) / sqrt_lambda)) * (1.0 - exp(-y / sqrt_lambda)) *
				(1.0 - exp((-y - 1.0) / sqrt_lambda)) * (3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) -
				2.0 * cos(4.0 * M_PI * x) * M_PI * (1.0 - exp(-x / sqrt_lambda)) / sqrt_lambda *
				exp((x - 1.0) / sqrt_lambda) * (1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				(3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) - 0.25 * (sin(4.0 * M_PI * x) + 2.0) *
				pow(sqrt_lambda, -2.0) * exp(-x / sqrt_lambda) * (1.0 - exp((x - 1.0) / sqrt_lambda)) * (1.0 - exp(-y / sqrt_lambda))
				* (1.0 - exp((-y - 1.0) / sqrt_lambda)) * (3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) - 0.5
				* (sin(4.0 * M_PI * x) + 2.0) * pow(sqrt_lambda, -2.0) * exp(-x / sqrt_lambda) * exp((x - 1.0) / sqrt_lambda)
				* (1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) * (3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) -
				0.25 * (sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) * pow(sqrt_lambda, -2.0) *
				exp((x - 1.0) / sqrt_lambda) * (1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				(3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) +
				(-0.25 * (sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) *
					(1.0 - exp((x - 1.0) / sqrt_lambda)) * exp(-y / sqrt_lambda) *
					(1.0 - exp((-y - 1.0) / sqrt_lambda)) * (3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) *
					pow(sqrt_lambda, -2.0) + 0.5 * (sin(4.0 * M_PI * x) + 2.0) *
					(1.0 - exp(-x / sqrt_lambda)) * (1.0 - exp((x - 1.0) / sqrt_lambda)) *
					exp(-y / sqrt_lambda) * exp((-y - 1.0) / sqrt_lambda) * (3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) *
					pow(sqrt_lambda, -2.0) - 0.5 * (sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) *
					(1.0 - exp((x - 1.0) / sqrt_lambda)) * exp(-y / sqrt_lambda) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
					exp((2.0 * y - 1.0) / sqrt_lambda / 2.0) * pow(sqrt_lambda, -2.0) - 0.25 *
					(sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) * (1.0 - exp((x - 1.0) / sqrt_lambda)) *
					(1.0 - exp(-y / sqrt_lambda)) * exp((-y - 1.0) / sqrt_lambda) * (3.0 - exp((2.0 * y - 1.0) / sqrt_lambda / 2.0)) *
					pow(sqrt_lambda, -2.0) - 0.5 * (sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) *
					(1.0 - exp((x - 1.0) / sqrt_lambda)) * (1.0 - exp(-y / sqrt_lambda)) * exp((-y - 1.0) / sqrt_lambda) *
					exp((2.0 * y - 1.0) / sqrt_lambda / 2.0) * pow(sqrt_lambda, -2.0) - 0.25 *
					(sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) * (1.0 - exp((x - 1.0) / sqrt_lambda)) *
					(1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) * exp((2.0 * y - 1.0) / sqrt_lambda / 2.0) *
					pow(sqrt_lambda, -2.0))) + 2 * (3 - exp((2 * y - 1) / (2 * sqrt_lambda))) * u_tmp; break;
			case 1:	return -lambda * (-4.0 * sin(4.0 * M_PI * x) * M_PI * M_PI * (1.0 - exp(-x / sqrt_lambda)) *
				(1.0 - exp((x - 1.0) / sqrt_lambda)) * (1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				(1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) + 0.200e1 * cos(4.0 * M_PI * x) *
				M_PI / sqrt_lambda * exp(-x / sqrt_lambda) * (1.0 - exp((x - 1.0) / sqrt_lambda)) *
				(1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				(1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) - 0.200e1 * cos(4.0 * M_PI * x) *
				M_PI * (1.0 - exp(-x / sqrt_lambda)) / sqrt_lambda * exp((x - 1.0) / sqrt_lambda) *
				(1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				(1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) - 0.25 * (sin(4.0 * M_PI * x) + 2.0) *
				pow(sqrt_lambda, -2.0) * exp(-x / sqrt_lambda) * (1.0 - exp((x - 1.0) / sqrt_lambda)) *
				(1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				(1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) - 0.5 * (sin(4.0 * M_PI * x) + 2.0) *
				pow(sqrt_lambda, -2.0) * exp(-x / sqrt_lambda) * exp((x - 1.0) / sqrt_lambda) *
				(1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				(1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) - 0.25 * (sin(4.0 * M_PI * x) + 2.0) *
				(1.0 - exp(-x / sqrt_lambda)) * pow(sqrt_lambda, -2.0) * exp((x - 1.0) / sqrt_lambda) *
				(1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				(1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) +
				-0.25 * (sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) *
				(1.0 - exp((x - 1.0) / sqrt_lambda)) * pow(sqrt_lambda, -2.0) * exp(-y / sqrt_lambda) *
				(1.0 - exp((-y - 1.0) / sqrt_lambda)) * (1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) +
				0.5 * (sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) *
				(1.0 - exp((x - 1.0) / sqrt_lambda)) * pow(sqrt_lambda, -2.0) * exp(-y / sqrt_lambda) *
				exp((-y - 1.0) / sqrt_lambda) * (1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) - 0.5 *
				(sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) *
				(1.0 - exp((x - 1.0) / sqrt_lambda)) * pow(sqrt_lambda, -2.0) * exp(-y / sqrt_lambda) *
				(1.0 - exp((-y - 1.0) / sqrt_lambda)) * exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0) - 0.25 *
				(sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) *
				(1.0 - exp((x - 1.0) / sqrt_lambda)) * (1.0 - exp(-y / sqrt_lambda)) *
				pow(sqrt_lambda, -2.0) * exp((-y - 1.0) / sqrt_lambda) *
				(1.0 + exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) - 0.5 *
				(sin(4.0 * M_PI * x) + 2.0) * (1.0 - exp(-x / sqrt_lambda)) *
				(1.0 - exp((x - 1.0) / sqrt_lambda)) * (1.0 - exp(-y / sqrt_lambda)) *
				pow(sqrt_lambda, -2.0) * exp((-y - 1.0) / sqrt_lambda) *
				exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0) + 0.25 * (sin(4.0 * M_PI * x) + 2.0) *
				(1.0 - exp(-x / sqrt_lambda)) * (1.0 - exp((x - 1.0) / sqrt_lambda)) *
				(1.0 - exp(-y / sqrt_lambda)) * (1.0 - exp((-y - 1.0) / sqrt_lambda)) *
				pow(sqrt_lambda, -2.0) * exp((1.0 - 2.0 * y) / sqrt_lambda / 2.0)) +
				(1 + exp((1 - 2 * y) / (2 * sqrt_lambda))) * u_tmp; break;
			}
		}

		if (Testing_parameters::test == 6)
			switch (area_number)
			{
			case 0: return 1.0; break;
			case 1:	return 1.0; break;
			}

		if (Testing_parameters::test == 7)
			switch (area_number)
			{
			case 0: return 0.0; break;
			}

		if (Testing_parameters::test == 8)
			switch (area_number)
			{
			case 0: return 0.0; break;
			case 1:	return 0.0; break;
			}

		return 1.0;
	}

	double Parameters::calculate_u_analytic(int area_number, double x, double y)
	{
		if (Testing_parameters::test == 1)
			switch (area_number)
			{
			case 0: return x + y; break;
			case 1:	return x + y; break;
			}

		if (Testing_parameters::test == 2)
			switch (area_number)
			{
			case 0: return x * x + y * y; break;
			case 1:	return x * x + y * y; break;
			}

		if (Testing_parameters::test == 3)
			switch (area_number)
			{
			case 0: return x * (1 - x) * y * (1 - y); break;
			case 1:	return x * (1 - x) * y * (1 - y); break;
			}

		if (Testing_parameters::test == 4)
			switch (area_number)
			{
			case 0: return x * x * y * y; break;
			case 1:	return x * x * y * y; break;
			}

		if (Testing_parameters::test == 5)
		{
			double lambda = calculate_lambda(area_number);
			double sqrt_lambda = sqrt(lambda);
			double u_tmp = 0.25 * (sin(4 * M_PI * x) + 2) * (1 - exp(-x / sqrt_lambda)) *
				(1 - exp((x - 1) / sqrt_lambda)) * (1 - exp(-y / sqrt_lambda))
				* (1 - exp((-y - 1) / sqrt_lambda));

			switch (area_number)
			{
			case 0: return (3 - exp((2 * y - 1) / (2 * sqrt_lambda))) * u_tmp; break;
			case 1:	return (1 + exp((1 - 2 * y) / (2 * sqrt_lambda))) * u_tmp; break;
			}
		}

		if (Testing_parameters::test == 6)
			switch (area_number)
			{
			case 0: return 0; break;
			case 1:	return 0; break;
			}

		if (Testing_parameters::test == 7)
			switch (area_number)
			{
			case 0: return 0; break;
			case 1:	return 0; break;
			}

		if (Testing_parameters::test == 8)
			switch (area_number)
			{
			case 0: return 0; break;
			case 1:	return 0; break;
			}

		return 1.0;
	}

	double Parameters::calculate_u_analytic2(partition::Partition p, double x, double y)
	{
		int el_num = p.search_element(x, y);
		return calculate_u_analytic(p.elements[el_num].number_of_area, x, y);
	}

	double Parameters::calculate_lambda(int area_number)
	{
		if (Testing_parameters::test == 1 || Testing_parameters::test == 2 || Testing_parameters::test == 3)
			switch (area_number)
			{
			case 0: return  1; break;
			case 1:	return  1; break;
			}

		if (Testing_parameters::test == 4)
			switch (area_number)
			{
			case 0: return 1; break;
			case 1:	return 1; break;
			}

		if (Testing_parameters::test == 5)
			switch (area_number)
			{
			case 0: return 1e-6; break;
			case 1:	return 1e-6; break;
			}
		if (Testing_parameters::test == 6)
			switch (area_number)
			{
			case 0:
			{
				//switch (Testing_parameters::test2)
				//{
				//	//case 1: return 1; break;//1-1å-9
				//case 1: return 1; break;//1-1å-9
				//case 2: return 1e-3; break;//1-1å-9
				//case 3: return 1e-5; break;//1-1å-9
				//case 4: return 1e-7; break;//1-1å-9
				//default: return 1; break;
				//}
				return 1e-3; break;//1-1å-9
			} break;
			}
		if (Testing_parameters::test == 7)
			switch (area_number)
			{
			case 0:
			{
				//switch (Testing_parameters::test2)
				//{
				//case 1: return 1; break;//1-1å-9
				//case 2: return 1e-3; break;//1-1å-9
				//case 3: return 1e-5; break;//1-1å-9
				//case 4: return 1e-7; break;//1-1å-9
				//default: return 1; break;
				//}
				return 1e-3; break;//1-1å-9
			} break;
			}

		if (Testing_parameters::test == 8)
			switch (area_number)
			{
			case 0:
			{
				//switch (Testing_parameters::test2)
				//{
				//case 1: return 1e-3; break;//1e-3-1å-6
				//case 2: return 1e-4; break;//1e-3-1å-6
				//case 3: return 1e-5; break;//1e-3-1å-6
				//case 4: return 1e-6; break;//1e-3-1å-6
				//default: return 1e-3; break;
				//}
				return 1e-3; break;//1-1å-9
			} break;

			case 1:
			{
				//switch (test2)
				//{
				//case 1: return 1e-3; break;//1e-3-1å-6
				//case 2: return 1e-4; break;//1e-3-1å-6
				//case 3: return 1e-5; break;//1e-3-1å-6
				//case 4: return 1e-6; break;//1e-3-1å-6
				//default: return 1e-3; break;
				//}
				return 1e-6; break;//1-1å-9
			} break;
			}

		return 1.0;
	}

	double Parameters::calculate_gamma(int area_number)
	{
		if (Testing_parameters::test == 1 || Testing_parameters::test == 2 || Testing_parameters::test == 3 || Testing_parameters::test == 4)
			switch (area_number)
			{
			case 0: return 0.0; break;
			case 1:	return 0.0; break;
			}
		if (Testing_parameters::test == 5)
			switch (area_number)
			{
			case 0: return 2.0; break;
			case 1:	return 1.0; break;
			}
		if (Testing_parameters::test == 6)
			switch (area_number)
			{
			case 0: return 10.0; break; //1-10
			case 1:	return 10.0; break;
			}

		if (Testing_parameters::test == 7)
			switch (area_number)
			{
			case 0: return 10; break;//1-10
			case 1:	return 10; break;
			}
		if (Testing_parameters::test == 8)
			switch (area_number)
			{
			case 0: return 0; break;
			case 1:	return 0; break;
			}

		return 1.0;
	}
	Point Parameters::calculate_a(int area_number)
	{
		Point a;
		if (Testing_parameters::test != 4 && Testing_parameters::test != 8) 
		{ 
			a.x = 0; 
			a.y = 0;
		}
		else 
		{ 
			if (Testing_parameters::test == 8) 
			{ 
				a.x = 1;
				a.y = 1;
			}
			else 
			{ 
				a.x = 1; 
				a.y = 1; 
			}
		}
		return a;
	}
}