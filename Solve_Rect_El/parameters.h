#pragma once
#define _USE_MATH_DEFINES
#include "partition.h"
namespace parameters
{

	class Parameters
	{
		public:
		
		//������ ������� �������
		static double Ug(int formula_number, double x, double y);
		//������� ������ �����
		static double calculate_f(int area_number, double x, double y);

		//������������� �������
		static double calculate_u_analytic(int area_number, double x, double y);
		static double calculate_u_analytic2(partition::Partition p, double x, double y);

		//������������ ������
		static double calculate_lambda(int area_number);
		static double calculate_gamma(int area_number);
		static point::Point calculate_a(int area_number);

	};

}	
