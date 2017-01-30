#pragma once
#include "partition.h"
#include "parameters.h"
#include "integration.h"
#include "basis.h"
#include "matrix.h"

namespace boundaries
{
	enum Formulation
	{BaumannOden, NIPG, IP};
	//численные потоки по внутренним границам
	class InternalBoundaries : virtual public partition::Partition,
							   virtual public parameters::Parameters,
							   virtual public integration::Gauss_integration,
							   virtual public basis::Basis
	{
		double mu; //коэффициенты стабилизации
		Formulation form;

		void calculate_E_horizontal(int element_number1, int element_number2, matrix::Matrix& A);
		void add_E_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector<std::vector<double>> &E);
		void calculate_E_vertical(int element_number1, int element_number2, matrix::Matrix& A);

	protected:

		void initialize(Formulation formulation);
		void calculate_internal_boundaries(int element_number, matrix::Matrix& A);
	};

	//численные потоки по внешним границам
	class OuterBoundaries : virtual public partition::Partition,
							virtual public parameters::Parameters,
							virtual public integration::Gauss_integration,
							virtual public basis::Basis
	{

		double mu; //коэффициенты стабилизации	
		Formulation form;

		void calculate_ES_out_left(int element_number, matrix::Matrix& A);
		void calculate_ES_out_right(int element_number, matrix::Matrix& A);
		void calculate_ES_out_low(int element_number, matrix::Matrix& A);
		void calculate_ES_out_up(int element_number, matrix::Matrix& A);

	protected:

		void initialize(Formulation formulation);
		void calculate_outer_boundaries(int element_number, matrix::Matrix& A);
	};

}