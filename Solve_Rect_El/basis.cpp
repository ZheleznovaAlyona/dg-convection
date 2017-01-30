#include "basis.h"
#include <array>
#include "myfunctions.h"

using namespace std;
using namespace partition;

namespace basis
{
	

	void Basis::initialize()
	{

		#pragma region непрерывный базис
		//n_func = 9;


		//array <function<double(double)>, n_func_1D> phi_;
		//phi_[0] = [](double ksi) { return 2 * (ksi - 0.5) * (ksi - 1); };
		//phi_[1] = [](double ksi) { return -4 * ksi * (ksi - 1); };
		//phi_[2] = [](double ksi) { return 2 * ksi * (ksi - 0.5); };

		//array <function<double(double)>, n_func_1D> dphi_ksi;
		//dphi_ksi[0] = [](double ksi) { return 4 * ksi - 3; };
		//dphi_ksi[1] = [](double ksi) { return  -8 * ksi + 4; };
		//dphi_ksi[2] = [](double ksi) { return 4 * ksi - 1; };

		//phi[0] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[0](etta); };
		//phi[1] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[0](etta); };
		//phi[2] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[0](etta); };
		//phi[3] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[1](etta); };
		//phi[4] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[1](etta); };
		//phi[5] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[1](etta); };
		//phi[6] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[2](etta); };
		//phi[7] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[2](etta); };
		//phi[8] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[2](etta); };

		//dphiksi[0] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[0](etta); };
		//dphiksi[1] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[0](etta); };
		//dphiksi[2] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[0](etta); };
		//dphiksi[3] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[1](etta); };
		//dphiksi[4] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[1](etta); };
		//dphiksi[5] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[1](etta); };
		//dphiksi[6] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[2](etta); };
		//dphiksi[7] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[2](etta); };
		//dphiksi[8] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[2](etta); };

		//dphietta[0] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[0](etta); };
		//dphietta[1] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[0](etta); };
		//dphietta[2] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[0](etta); };
		//dphietta[3] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[1](etta); };
		//dphietta[4] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[1](etta); };
		//dphietta[5] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[1](etta); };
		//dphietta[6] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[2](etta); };
		//dphietta[7] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[2](etta); };
		//dphietta[8] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[2](etta); };
#pragma endregion


		#pragma region разрывный базис
		array <function<double(double)>, n_func_1D> phi_;
		phi_[0] = [](double ksi) { return 0.5 * (ksi - 1) * (1 - 3 * ksi); };
		phi_[1] = [](double ksi) { return -2 * pow_i(2, ksi) + 2 * ksi; };
		phi_[2] = [](double ksi) { return 0.5 * (3 * pow_i(2, ksi) - 2 * ksi); };

		array <function<double(double)>, n_func_1D> dphi_ksi;
		dphi_ksi[0] = [](double ksi) { return -3 * ksi + 2; };
		dphi_ksi[1] = [](double ksi) { return  -4 * ksi + 2; };
		dphi_ksi[2] = [](double ksi) { return 3 * ksi - 1; };

		phi[0] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[0](etta); };
		phi[1] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[0](etta); };
		phi[2] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[0](etta); };
		phi[3] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[1](etta); };
		phi[4] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[1](etta); };
		phi[5] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[1](etta); };
		phi[6] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[2](etta); };
		phi[7] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[2](etta); };
		phi[8] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[2](etta); };

		dphiksi[0] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[0](etta); };
		dphiksi[1] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[0](etta); };
		dphiksi[2] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[0](etta); };
		dphiksi[3] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[1](etta); };
		dphiksi[4] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[1](etta); };
		dphiksi[5] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[1](etta); };
		dphiksi[6] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[2](etta); };
		dphiksi[7] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[2](etta); };
		dphiksi[8] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[2](etta); };

		dphietta[0] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[0](etta); };
		dphietta[1] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[0](etta); };
		dphietta[2] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[0](etta); };
		dphietta[3] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[1](etta); };
		dphietta[4] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[1](etta); };
		dphietta[5] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[1](etta); };
		dphietta[6] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[2](etta); };
		dphietta[7] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[2](etta); };
		dphietta[8] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[2](etta); };

		#pragma endregion

		#pragma region линейный базис
		//array <function<double(double)>, n_func_1D> phi_;
		//phi_[0] = [](double ksi) { return 1 - ksi; };
		//phi_[1] = [](double ksi) { return ksi; };

		//array <function<double(double)>, n_func_1D> dphi_ksi;
		//dphi_ksi[0] = [](double ksi) { return -1; };
		//dphi_ksi[1] = [](double ksi) { return  1; };

		//phi[0] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[0](etta); };
		//phi[1] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[0](etta); };
		//phi[2] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[1](etta); };
		//phi[3] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[1](etta); };

		//dphiksi[0] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[0](etta); };
		//dphiksi[1] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[0](etta); };
		//dphiksi[2] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[1](etta); };
		//dphiksi[3] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[1](etta); };

		//dphietta[0] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[0](etta); };
		//dphietta[1] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[0](etta); };
		//dphietta[2] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[1](etta); };
		//dphietta[3] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[1](etta); };

		#pragma endregion

		#pragma region базис в другом виде
		//phi[0] = [](double ksi, double etta)
		//{ return 0.5 * (ksi - 1) * (1 - 3 * ksi) * 0.5 * (etta - 1) * (1 - 3 * etta); };
		//phi[1] = [](double ksi, double etta)
		//{ return (-2 * ksi * ksi + 2 * ksi) * 0.5 * (etta - 1) * (1 - 3 * etta); };
		//phi[2] = [](double ksi, double etta)
		//{ return 0.5 * (3 * ksi * ksi - 2 * ksi) * 0.5 * (etta - 1) * (1 - 3 * etta); };
		//phi[3] = [](double ksi, double etta)
		//{ return 0.5 * (ksi - 1) * (1 - 3 * ksi) * (-2 * etta * etta + 2 * etta); };
		//phi[4] = [](double ksi, double etta)
		//{ return (-2 * ksi * ksi + 2 * ksi) * (-2 * etta * etta + 2 * etta); };
		//phi[5] = [](double ksi, double etta)
		//{ return 0.5 * (3 * ksi * ksi - 2 * ksi) * (-2 * etta * etta + 2 * etta); };
		//phi[6] = [](double ksi, double etta)
		//{ return 0.5 * (ksi - 1) * (1 - 3 * ksi) * 0.5 * (3 * etta * etta - 2 * etta); };
		//phi[7] = [](double ksi, double etta)
		//{ return (-2 * ksi * ksi + 2 * ksi) * 0.5 * (3 * etta * etta - 2 * etta); };
		//phi[8] = [](double ksi, double etta)
		//{ return 0.5 * (3 * ksi * ksi - 2 * ksi) * 0.5 * (3 * etta * etta - 2 * etta); };

		//dphiksi[0] = [](double ksi, double etta)
		//{ return (0.5 - 1.5 * ksi) * (0.5 * etta - 0.5) * (1 - 3 * etta)
		//	- (1.5 * ksi - 1.5) * (0.5 * etta - 0.5) * (1 - 3 * etta); };
		//dphiksi[1] = [](double ksi, double etta)
		//{ return (-4 * ksi + 2) * (0.5 * etta - 0.5) * (1 - 3 * etta); };
		//dphiksi[2] = [](double ksi, double etta)
		//{ return (3 * ksi - 1) * (0.5 * etta - 0.5) * (1 - 3 * etta); };
		//dphiksi[3] = [](double ksi, double etta)
		//{ return (0.5 - 1.5 * ksi) * (-2 * etta * etta + 2 * etta)
		//	- (1.5 * ksi - 1.5) * (-2 * etta * etta + 2 * etta); };
		//dphiksi[4] = [](double ksi, double etta)
		//{ return (-4 * ksi + 2) * (-2 * etta * etta + 2 * etta); };
		//dphiksi[5] = [](double ksi, double etta)
		//{ return (3 * ksi - 1) * (-2 * etta * etta + 2 * etta); };
		//dphiksi[6] = [](double ksi, double etta)
		//{ return (0.5 - 1.5 * ksi) * (1.5 * etta * etta - etta)
		//	- (1.5 * ksi - 1.5) * (1.5 * etta * etta - etta); };
		//dphiksi[7] = [](double ksi, double etta)
		//{ return (-4 * ksi + 2) * (1.5 * etta * etta - etta); };
		//dphiksi[8] = [](double ksi, double etta)
		//{ return (3 * ksi - 1) * (1.5 * etta * etta - etta); };

		//dphietta[0] = [](double ksi, double etta)
		//{ return (0.25 * ksi - 0.25) * (1 - 3 * ksi) * (1 - 3 * etta)
		//	- (1.5 * ksi - 1.5) * (1 - 3 * ksi) * (0.5 * etta - 0.5); };
		//dphietta[1] = [](double ksi, double etta)
		//{ return (-ksi * ksi + ksi) * (1 - 3 * etta)
		//	- (-6 * ksi * ksi + 6 * ksi) * (0.5 * etta - 0.5); };
		//dphietta[2] = [](double ksi, double etta)
		//{ return (0.75 * ksi * ksi - 0.5 * ksi) * (1 - 3 * etta)
		//	- (4.5 * ksi * ksi - 3 * ksi) * (0.5 * etta - 0.5); };
		//dphietta[3] = [](double ksi, double etta)
		//{ return (0.5 * ksi - 0.5) * (1 - 3 * ksi) * (-4 * etta + 2); };
		//dphietta[4] = [](double ksi, double etta)
		//{ return (-2 * ksi * ksi + 2 * ksi) * (-4 * etta + 2); };
		//dphietta[5] = [](double ksi, double etta)
		//{ return (1.5 * ksi * ksi - ksi) * (-4 * etta + 2); };
		//dphietta[6] = [](double ksi, double etta)
		//{ return (0.5 * ksi - 0.5) * (1 - 3 * ksi) * (3 * etta - 1); };
		//dphietta[7] = [](double ksi, double etta)
		//{ return (-2 * ksi * ksi + 2 * ksi) * (3 * etta - 1); };
		//dphietta[8] = [](double ksi, double etta)
		//{ return (1.5 * ksi * ksi - ksi) * (3 * etta - 1); };
		#pragma endregion


	}

	double Basis::phi_i(int i, double x, double y, int element_number, Partition& p)
	{
		double x_left = p.nodes[p.elements[element_number].nodes[0]].x;
		double x_right = p.nodes[p.elements[element_number].nodes[1]].x;
		double y_low = p.nodes[p.elements[element_number].nodes[0]].y;
		double y_up = p.nodes[p.elements[element_number].nodes[3]].y;
		double hx = x_right - x_left, hy = y_up - y_low;
		double ksi = (x - x_left) / hx, etta = (y - y_low) / hy;
		return phi[i](ksi, etta);
	}
}