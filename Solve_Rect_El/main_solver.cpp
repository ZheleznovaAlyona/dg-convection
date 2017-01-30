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
	MainSolver::MainSolver() {}

	MainSolver::MainSolver(std::ifstream& grid_f_in,
		std::ifstream& elements_f_in,
		std::string log_f,
		std::ifstream& boundary1,
		std::ifstream& boundary2,
		std::ifstream& boundary3,
		Formulation formulation)
	{
		initialize(grid_f_in,
			elements_f_in,
			log_f,
			boundary1,
			boundary2,
			boundary3,
			formulation);
	}

	MainSolver::~MainSolver() {}

	void MainSolver::initialize(ifstream& grid_f_in,
								ifstream& elements_f_in,
								string log_f,
								ifstream& boundary1,
								ifstream& boundary2,
								ifstream& boundary3,
								Formulation formulation)
	{
		boundary1 >> boundaries1;
		boundary2 >> boundaries2;
		boundary3 >> boundaries3;
		BoundaryConditionsSupport::form_ = formulation;

		Partition::input(grid_f_in, elements_f_in);
		Basis::initialize();
		Gauss_integration::initialize();
		int slae_size = 0;
		for (unsigned int i = 0; i < elements.size(); i++)
			slae_size += elements[i].ndof;

		InternalBoundaries::initialize(formulation);
		OuterBoundaries::initialize(formulation);
		//BoundaryConditionsSupport::initialize_penalty_parameters();

		logger.open(log_f);

		int tmp = Partition::count_unzero_matrix_elements(slae_size);
		my_slae.initialize(slae_size, tmp);

		U_numerical.initialize(nodes.size());
		q_prev.initialize(slae_size);
	}

	void MainSolver::reinitialize()
	{
		U_numerical.make_zero();
		my_slae.reinitialize();
	}

	void MainSolver::build_slae(MyVector q_calc)
	{
		logger.send_message_build_slae();
		int size = Partition::elements.size();

		//локаьные матрицы и вектор правой части
		for(int el_i = 0; el_i < size; el_i++)
			calculate_locals(el_i, q_calc);

		//матрицы межэлементных границ
		for(int el_i = 0; el_i < size; el_i++)
			InternalBoundaries::calculate_internal_boundaries(el_i, my_slae.A);

		//расчёт матриц границ области
			for(int el_i = 0; el_i < size; el_i++)
				OuterBoundaries::calculate_outer_boundaries(el_i, my_slae.A);

		//учёт первых краевых условий
		BoundaryConditionsSupport::calculate_all_boundaries1(my_slae.b);
	}

	void MainSolver::linear(ofstream& solution_f_out,
							ofstream& info_f_out,
							Solver& s)
	{
		//bool one_research = true;
		double normL2u;
		MyVector sol_0, U_an;
		U_an.initialize(nodes.size());
		sol_0.initialize(my_slae.n);

		my_slae.A.create_portret(static_cast<Partition>(*this), logger);

		//vector<Point> nd;
		//nd.reserve(my_slae.n);
		//for (size_t i = 0; i < elements.size(); i++)
		//{
		//	Point nds[9];
		//	nds[0] = nodes[elements[i].nodes[0]];
		//	nds[2] = nodes[elements[i].nodes[1]];
		//	nds[6] = nodes[elements[i].nodes[2]];
		//	nds[8] = nodes[elements[i].nodes[3]];

		//	nds[1] = (nds[0] + nds[2]) * 0.5;
		//	nds[3] = (nds[0] + nds[6]) * 0.5;
		//	nds[5] = (nds[2] + nds[8]) * 0.5;
		//	nds[7] = (nds[6] + nds[8]) * 0.5;

		//	nds[4] = (nds[1] + nds[7]) * 0.5;
		//	for (size_t j = 0; j < 9; j++)
		//	{
		//		nd.push_back(nds[j]);
		//	}

		//}

		//for (size_t i = 0; i < my_slae.n; i++)
		//{
		//	sol_0[i] = calculate_u_analytic2(static_cast<Partition>(*this), nd[i].x, nd[i].y);
		//}

		build_slae(sol_0);



		//ofstream ff("matrix.txt");
		//for (size_t i = 0; i < my_slae.n; i++)
		//{
		//	for (size_t j = 0; j < my_slae.n; j++)
		//	{
		//		ff << my_slae.A.mtr[i][j] << "\t";
		//	}
		//	ff << endl;
		//}

		//ofstream ff2("matrixggu.txt");
		//for (size_t i = 0; i < my_slae.A.size; i++)
		//{
		//	ff2 << my_slae.A.ggu[i] << endl;
		//}

		//ofstream ff3("matrixggl.txt");
		//for (size_t i = 0; i < my_slae.A.size; i++)
		//{
		//	ff3 << my_slae.A.ggl[i] << endl;
		//}

		//ofstream ff4("matrixdi.txt");
		//for (size_t i = 0; i < my_slae.n; i++)
		//{
		//	ff4 << my_slae.A.di[i] << endl;
		//}

		logger.send_message_solution();
		q_prev = s.solve(sol_0, normL2u, my_slae, logger);

		get_vector_solution_in_nodes_u(q_prev, U_numerical);
		get_vector_in_nodes_u_an(U_an);
		

		normL2u = diff_normL2_u(q_prev); 

		logger.output(solution_f_out, info_f_out, normL2u, U_numerical, static_cast<Partition>(*this), U_an);
		string name = "out/solution.dat";
		out_tecplot(name, q_prev);
	}

void MainSolver::solve(std::ofstream & solution_f_out, std::ofstream & info_f_out)
{
	MyVector q(my_slae.n);
	Solver *s;

	switch (Testing_parameters::solver)
	{
	case 1:
	{
		s = new BiCGStab();
	}
	break;
	case 2:
	{
		s = new GMRES();
	}
	break;
	case 3:
	{
		s = new BCGandGMRESSolver();
	}
	break;
	case 4:
	{
		s = new BCG();
	}
	break;
	default:
	{
		s = new GMRES();
	}
	};

	if(Testing_parameters::use_LU) my_slae.A.LU();

	s->s_parameters.initialize("in/solver.json");

	linear(solution_f_out, info_f_out, *s);
}

double MainSolver::diff_normL2_u(MyVector q_solution)
{
	double diff_local, diff;
	int size = elements.size();
	double u, function;
	double x0, y0, hx, hy;
	double jacobian;

	diff = 0;
	for(int i = 0; i < size; i++)
	{
		x0 = nodes[elements[i].nodes[0]].x;
		y0 = nodes[elements[i].nodes[0]].y;
		hx = get_hx(i);
		hy = get_hy(i);
		jacobian = hx * hy / 4.0;

		diff_local = 0;
		for(int k = 0; k < 9; k++)
		{
			double p_ksi = 0.5 + 0.5 * gauss_points[0][k],
				   p_etta = 0.5 + 0.5 * gauss_points[1][k];
			//double p_x = hx * gauss_points[0][k] + x0;
			//double p_y = hy * gauss_points[1][k] + y0;
			double p_x = hx * p_ksi + x0;
			double p_y = hy * p_etta + y0;
			u = get_solution_in_point_u(p_x, p_y, i, q_solution);
			function = calculate_u_analytic(elements[i].number_of_area, p_x, p_y);
			function -= u;
			diff_local += gauss_weights[k] * function * function;
		}
		diff_local *= jacobian;

		diff += diff_local;
	}

	return sqrt(diff / size);
}

void MainSolver::out_tecplot(string & tecplot_filename, MyVector& q)
{
	static auto& part = static_cast<Partition>(*this);
	double min_x = nodes[0].x, max_x = nodes[nodes.size() - 1].x, step_x = 0.01;
	double min_y = nodes[0].y, max_y = nodes[nodes.size() - 1].y, step_y = 0.01;


	size_t nx = (size_t)((max_x - min_x) / step_x);
	size_t ny = (size_t)((max_y - min_y) / step_y);

	printf("\nWriting to Tecplot ...\n");

	ofstream tecplot_file(tecplot_filename.c_str());
	ofstream tecplot_file2("out/an.dat");

	tecplot_file << "TITLE = \"Problem 2\"\n";
	tecplot_file << "VARIABLES = \"x\", \"y\", \"U\"\n";
	tecplot_file << "ZONE I= " << nx + 1 << ", J= " << ny + 1 << ", F=POINT\n";

	tecplot_file2 << "TITLE = \"Problem 2\"\n";
	tecplot_file2 << "VARIABLES = \"x\", \"y\", \"U\"\n";
	tecplot_file2 << "ZONE I= " << nx + 1 << ", J= " << ny + 1 << ", F=POINT\n";

	tecplot_file.precision(17);
	tecplot_file << scientific;

	tecplot_file2.precision(17);
	tecplot_file2 << scientific;

	for (size_t j = 0; j <= ny; j++)
	{
		double y = min_y + step_y * (double)j;
		for (size_t k = 0; k <= nx; k++)
		{
			double x = min_x + step_x * (double)k;

			Point p = Point(x, y);
			double U_sol = get_solution_in_point2_u(x, y, q);
			tecplot_file << x << " " << y << " " << U_sol << "\n";
			tecplot_file2 << x << " " << y << " " << calculate_u_analytic2(part, x, y) << "\n";
		}
	}

	tecplot_file << "\n";
	tecplot_file.flush();
	tecplot_file.close();

	tecplot_file2 << "\n";
	tecplot_file2.flush();
	tecplot_file2.close();
}

}
