/* Copyright (C) 2018 Bho Matthiesen
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
#include "HKSolver.h"

class HK_GEE : public HKSolver<>
{
	public:
		const double mu = 4.0;
		const double Pc = 1.0;

		template <typename... Ts>
		HK_GEE(Ts&&... args)
			: HKSolver<>(std::forward<Ts>(args)...)
		{
			typename HK::vtypeC linC = {};
			std::fill_n(linC.begin(), HKDIM, 1.0);
			
			typename HK::vtypeC denomLinC = {};
			std::fill_n(denomLinC.begin()+HKDIM, HKDIM, mu);

			typename HK::vtype denomLinNC = {};
			std::fill_n(denomLinNC.begin(), HKDIM, mu);

			s.setObjective(linC, {}, denomLinC, denomLinNC, 0.0, Pc);
			s.setPrecision(1e-3, 1e-5);
		}

		std::string getName() const override
		{
			return "hkEE";
		}
};

int main(int argc, char *argv[])
#ifndef DEBUG
	try
#endif
{
	Sim<HK_GEE> s(argc, argv);
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
