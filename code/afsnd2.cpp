/* Copyright (C) 2018-2019 Bho Matthiesen
 * 
 * This program is used in the article:
 * 
 * Bho Matthiesen and Eduard Jorswieck, "Efficient Global Optimal Resource
 * Allocation in Non-Orthogonal Interference Networks," submitted to IEEE
 * Transactions on Signal Processing.
 * 
 * 
 * License:
 * This program is licensed under the GPLv2 license. If you in any way use this
 * code for research that results in publications, please cite our original
 * article listed above.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. */

// SND
/* afsnd code starts here */
#include <iostream>
#include <stdexcept>

#include "gurobi_c++.h"

#include "Sim.h"
#include "SNDSolver.h"

int main(int argc, char *argv[])
#ifndef DEBUG
	try
#endif
{
	Sim<SNDSolver<>> s(argc, argv);
	s.run();

	return 0;
}
#ifndef DEBUG
catch (std::exception &e)
{
	std::cerr << "Error: " << e.what() << std::endl;
	return 1;
}
catch (GRBException &e)
{
	std::cerr << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")" << std::endl;
	return 1;
}
catch (...)
{
	std::cerr << "Error: Unknown exception" << std::endl;
	throw;
	return 1;
}
#endif
