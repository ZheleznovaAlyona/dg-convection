#include <fstream>
#include "main_solver.h"
#include "testing_parameters.h"
#include "boundaries.h"

using namespace std;
using namespace mainsolver;
using namespace solver;
using namespace testingparameters;
using namespace boundaries;

void main()
{
	ifstream l1_in("in/l1.dt"), l2_in("in/l2.dt"), l3_in("in/l3.dt"), grid_in("in/grid.dt"), elements_in("in/elements.dt");
	ofstream sol_out("out/solution.txt"), info_out("out/info.txt");
	Testing_parameters::initialize("in/testing_parameters.json");

	MainSolver problem_solver(grid_in, elements_in, "out/log.txt", l1_in, l2_in, l3_in, IP);	
	problem_solver.solve(sol_out, info_out);

	l1_in.close();
	l2_in.close();
	l3_in.close();
	grid_in.close();
	elements_in.close();

	system("pause");
}