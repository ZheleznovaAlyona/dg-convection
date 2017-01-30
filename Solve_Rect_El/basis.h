#pragma once
#include <functional>
#include "partition.h"

namespace basis
{
	const int n_func = 9;
	const int n_func_1D = 3;

	//const int n_func = 4;
	//const int n_func_1D = 2;

	class Basis
	{

	protected:

		//int n_func; //число базисных функций

		std::function<double(double, double)> phi[n_func]; //указатели на функции вычисления базисных функций в точке
		std::function<double(double, double)> dphiksi[n_func]; //указатели на функции вычисления d/dksi базисных функций в точке
		std::function<double(double, double)> dphietta[n_func]; //указатели на функции вычисления d/detta базисных функций в точке

		void initialize();
		double phi_i(int i, double x, double y, int element_number, partition::Partition& p);
	};
}