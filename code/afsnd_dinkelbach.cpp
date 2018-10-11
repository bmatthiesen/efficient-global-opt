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
#include "SNDSolver.h"

#ifndef CASE
#define CASE 7
#endif

#define useTighterBound false

using _SND = SND<CASE>;

class Dinkelbach
{
	static_assert(CASE <= 7, "Out of range");

	public:
		using H5 = HDF<SNDDIM,SNDDIM>;

		const double outerTol = 1e-3;
		const double mu = 4.0;
		const double Pc = 1.0;

		Dinkelbach(H5::WPPtr wp) : wp(wp) { }

		const std::string getName() const
		{
			return "dinkelbach";
		}

		constexpr size_t getDim()
		{
			return 1;
		}

		void enableBackup(const std::string&, const std::chrono::seconds&)
			{ std::cerr << "Backup not implemented" << std::endl; }

		void removeBackup()
			{ }

		const typename H5::Result solve(const size_t idx);

	private:
		H5::WPPtr wp;

		using clk = std::chrono::high_resolution_clock;

		double lambda = 0.0;
};


const typename Dinkelbach::H5::Result
Dinkelbach::solve(const size_t)
{
		unsigned long long iter = 0;
		unsigned long long innerIter = 0;
		double cbv = 0;
		std::chrono::time_point<clk> tic = clk::now();

		typename _SND::vtypeC linC;
		typename _SND::vtype linNC;

		linC.fill(1.0);
		lambda = 0.0;

		// dinkelbach
		std::unique_ptr<_SND> s;
		do
		{
			s = std::make_unique<_SND>();
			++iter;

			s->init(wp, useTighterBound);

			linNC.fill(-1.0*mu*lambda);
			s->setObjective(linC, linNC, {}, {}, -1.0 * lambda * Pc);

			// call solver
			try
			{
				s->optimize();
			}
			catch (GRBException &e)
			{
				std::stringstream ss;
				ss << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")";
				throw std::runtime_error(ss.str());
			}

			cbv = s->optval;
			double numer = std::accumulate(s->optC.begin(), s->optC.end(), 0.0);
			double denom = mu * std::accumulate(s->optNC.begin(), s->optNC.end(), 0.0) + Pc;
			lambda = numer / denom;

			innerIter += s->iter;
		} while (std::abs(cbv) > outerTol);

		// output
		std::chrono::time_point<clk> toc = clk::now();

		H5::Result res(s.get());
		res.objective = lambda;
		res.iter = iter;
		res.lastUpdate = innerIter;
		res.runtime = std::chrono::duration<double>(toc - tic).count();

		std::cout << std::endl << std::endl;
		std::cout << "lambda: " << lambda << std::endl;
		std::cout << "iter: " << iter << std::endl;
		std::cout << "total inner iter: " << innerIter << std::endl;
		std::cout << "Time: " << res.runtime << std::endl;

		return res;
}


int main(int argc, char *argv[])
#ifndef DEBUG
	try
#endif
{
	Sim<Dinkelbach> s(argc, argv);
	s.run();

	return 0;
}
#ifndef DEBUG
catch (std::exception &e)
{
	std::cerr << "Error: " << e.what() << std::endl;
	return 1;
}
catch (...)
{
	std::cerr << "Error: Unknown exception" << std::endl;
	throw;
	return 1;
}
#endif
