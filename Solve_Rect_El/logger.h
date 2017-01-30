#pragma once
#include <fstream>
#include <string>
#include "myvector.h"
#include "partition.h"
#include "parameters.h"

namespace logger
{
	class Logger
	{
		std::ofstream log_f;

		enum 
		{
			TO_FILE,
			TO_STDO
		} p;

	public:
		void open(std::string file_name);
		Logger();
		~Logger();

		void send_current_information(double r_norm, int iteration_number);
		void send_current_information_to_screen(double r_norm, int iteration_number);

		void send_message_build_slae();
		void send_message_create_portret();
		void send_message_solution();
		void send_message_U();

		void output(std::ofstream& solution_f_out, 
					std::ofstream& info_f_out, 
					double normL2u, 
					myvector::MyVector& U,
					partition::Partition& P,
					myvector::MyVector& U_an);

	};
}