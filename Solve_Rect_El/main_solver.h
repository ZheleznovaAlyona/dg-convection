#pragma once
#include <string>
#include <fstream>

#include "boundary_conditions.h"
#include "myvector.h"
#include "boundaries.h"
#include "logger.h"
#include "solver.h"
#include "SLAE.h"

namespace mainsolver
{
	class MainSolver : public boundary_conditions::BoundaryConditionsSupport, public boundaries::InternalBoundaries, public boundaries::OuterBoundaries
	{
		slae::SLAE my_slae;
		logger::Logger logger; //логгер для вывода информации о процессе решения СЛАУ

		myvector::MyVector U_numerical; //численное решение
		myvector::MyVector q_prev; //веса

		void initialize(std::ifstream& grid_f_in,
			std::ifstream& elements_f_in,
			std::string log_f,
			std::ifstream& boundary1,
			std::ifstream& boundary2,
			std::ifstream& boundary3,
			boundaries::Formulation formulation);

		void reinitialize();

		void build_slae(myvector::MyVector q_calc);

		void linear(std::ofstream& solution_f_out,
					std::ofstream& info_f_out,
					solver::Solver& s);
		double get_solution_in_point_u(double x, double y, int element_number, myvector::MyVector qi);
		double get_solution_in_point2_u(double x, double y, myvector::MyVector qi);
		void get_vector_solution_in_nodes_u(myvector::MyVector qi, myvector::MyVector &solution);
		void get_vector_in_nodes_u_an(myvector::MyVector &U_an);

		//локальные матрицы и векторы
		void calculate_locals(int element_number, myvector::MyVector q_calc);
		void calculate_G(int element_number);
		void calculate_C(int element_number, myvector::MyVector q_calc);
		void calculate_M(int element_number);
		void calculate_F(int element_number);

		double diff_normL2_u(myvector::MyVector q_solution);//погрешность решения в норме L2
		void out_tecplot(std::string & tecplot_filename, myvector::MyVector& q);

	public:
   
		MainSolver();

		MainSolver(std::ifstream& grid_f_in,
				   std::ifstream& elements_f_in,
				   std::string log_f,
				   std::ifstream& boundary1,
				   std::ifstream& boundary2,
				   std::ifstream& boundary3,
				   boundaries::Formulation formulation);
					 
		~MainSolver();
		void solve(std::ofstream& solution_f_out, std::ofstream& info_f_out);							   
	};
}