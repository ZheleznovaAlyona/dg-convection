#include "solver.h"
#include "myfunctions.h"

using namespace myvector;
using namespace logger;

namespace solver
{
	MyVector BCGandGMRESSolver::solve(MyVector& U_begin, double &normL2u, slae::SLAE& slae_in, Logger& my_logger)
	{
		MyVector q(slae_in.n);
			
		double eps2 = GMRES::s_parameters.epsilon;
		GMRES::s_parameters.epsilon = 1e-7;
		q = GMRES::solve(U_begin, normL2u, slae_in, my_logger);
		GMRES::s_parameters.epsilon = eps2;
		q = BCG::solve(q, normL2u, slae_in, my_logger);
		return q;
	}
}