#include "logger.h"
#include <iomanip>
#include <iostream>

using namespace std;
using namespace myvector;
using namespace partition;
using namespace parameters;

namespace logger
{
	void Logger::send_current_information(double r_norm, int iteration_number)
	{
		ostream* os;

		if(p == TO_FILE)
			os = &log_f;
		else
			os = &cout;

		*os << iteration_number << "     " << scientific << setprecision(20) << r_norm << endl;
	}

	void Logger::send_current_information_to_screen(double r_norm, int iteration_number)
	{
		cout << iteration_number << "\tr=" << scientific << setprecision(10) << r_norm << endl;
	}

	void Logger::send_message_build_slae()
	{
		cout << "Building SLAE..." << endl;
	}

	void Logger::send_message_create_portret()
	{
		cout << "Creating matrix-portret..." << endl;
	}

	void Logger::send_message_solution()
	{
		cout << "Getting solution..." << endl;
	}

	void Logger::send_message_U()
	{
		cout << "Getting U..." << endl;
	}

	Logger::Logger()
	{
		p = TO_STDO;
	}

	void Logger::open(string file_name)
	{
		log_f.open(file_name);
		p = TO_FILE;
	}

	Logger::~Logger()
	{
		if(p == TO_FILE)
			log_f.close();
	}


	void Logger::output(ofstream& solution_f_out, 
						ofstream& info_f_out, 
						double normL2u, 
						MyVector& U,
						Partition& P,
						MyVector& U_an)
	{		
		ofstream an("out/an.txt");
		solution_f_out << U << endl << endl << endl;

		info_f_out << "norm L2 u:|u*-u|=" << scientific << setprecision(4) << normL2u 
			<< endl;
		an << U_an << endl << endl << endl;
	}
}