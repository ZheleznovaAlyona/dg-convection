#include "SLAE.h"

namespace slae
{

	void SLAE::initialize(int size, int unzero_matrix_elements)
	{
		n = size;
		A.initialize(n, unzero_matrix_elements);
		//A.mtr.resize(n);
		//for (size_t i = 0; i < n; i++)
		//{
		//	A.mtr[i].resize(n);
		//	memset(&A.mtr[i][0], 0, sizeof(double) * n);
		//}
		b.initialize(n);
	}

	void SLAE::reinitialize()
	{
		A.reinitialize();
		b.make_zero();
	}

}
