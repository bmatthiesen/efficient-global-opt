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

#include <iostream>
#include <stdexcept>

#include "Sim.h"
#include "SNDSolver.h"

#ifndef CASE
#define CASE 8 // 7 for classical SND, 0 for classic TIN, -1 for all SND cases
#endif

template <size_t C = 8, bool useTighterBound = false>
class GEE : public SNDSolver<C,useTighterBound>
{
	public:
		const double mu = 4.0;
		const double Pc = 1.0;

		template <typename... Ts>
		GEE(Ts&&... args)
			: SNDSolver<C,useTighterBound>(std::forward<Ts>(args)...)
		{
			typename SND<0>::vtypeC linC;
			linC.fill(1.0);
			
			typename SND<0>::vtype denomLinNC;
			denomLinNC.fill(mu);

			this->setObjective(linC, {}, {}, denomLinNC, 0.0, Pc);
			this->setPrecision(1e-3, 1e-5);
		}

		std::string getName() const override
		{
			if constexpr(C == 7)
				return "afsndEE";
			else if constexpr(C == 8)
				return "afsnd2EE";
			else if constexpr(C == 0)
				return "aftinEE";
			else
				return std::string("af") + std::to_string(C) + std::string("EE");
		}
};


int main(int argc, char *argv[])
#ifndef DEBUG
	try
#endif
{
#ifdef MOSEK
	Sim<GEE<CASE,true>> s(argc, argv, BndSolvers::Mosek);
#else
	Sim<GEE<CASE>> s(argc, argv);
#endif
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
