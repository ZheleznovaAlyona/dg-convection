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

		//int n_func; //����� �������� �������

		std::function<double(double, double)> phi[n_func]; //��������� �� ������� ���������� �������� ������� � �����
		std::function<double(double, double)> dphiksi[n_func]; //��������� �� ������� ���������� d/dksi �������� ������� � �����
		std::function<double(double, double)> dphietta[n_func]; //��������� �� ������� ���������� d/detta �������� ������� � �����

		void initialize();
		double phi_i(int i, double x, double y, int element_number, partition::Partition& p);
	};
}