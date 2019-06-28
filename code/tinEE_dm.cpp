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

#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <string>

#include <iostream>

#include "MMP.h"
#include "HDF.h"
#include "Sim.h"

using std::array;

namespace mmp {
template <size_t Dim>
class TIN_Dinkelbach : public MMP<Dim>
{
	using typename BRB<Dim>::vtype;

	public:
		using H5 = HDF<Dim,1>;

		double w[Dim];
		double alpha[Dim];
		array<double, Dim> beta[Dim];
		double sigma[Dim];

		const double mu = 4.0;
		const double outerTol = 1e-4;
		const double Pc = 1.0;

		TIN_Dinkelbach(typename H5::WPPtr wp);

		typename H5::Result solve(const size_t = 0);

		const std::string getName() const
		{
			return "tin_dinkelbach";
		}

		constexpr size_t getDim()
		{
			return 1;
		}

		void enableBackup(const std::string&, const std::chrono::seconds&)
			{ std::cerr << "Backup not implemented" << std::endl; }

		void removeBackup()
			{ }

	private:
		typename H5::WPPtr wp;

		// avoid overriding
		using MMP<Dim>::obj;
		using MMP<Dim>::isFeasible;

		double obj(const vtype& x, const vtype& y) const override;
		bool isFeasible(const vtype&, const vtype&) const override { return true; };

		using clk = std::chrono::high_resolution_clock;

		double lambda;

		template <typename T>
		static double abssquare(const T& a)
			{ return a.re * a.re + a.im * a.im; }
};


template <size_t D>
TIN_Dinkelbach<D>::TIN_Dinkelbach(typename H5::WPPtr wp)
	: MMP<D>(), wp(wp)
{
	const unsigned short mx_q[] = {1,2,0};
	const unsigned short mx_l[] = {2,0,1};

	auto P0 = std::pow(10, wp->P0/10.0);
	decltype(wp->P) P;
	std::transform(wp->P.begin(), wp->P.end(), P.begin(), [](const double p) { return std::pow(10, p/10.0);});

	typename H5::vec_t habs;
	typename H5::vec_t gabs;

	std::transform(wp->h.begin(), wp->h.end(), habs.begin(), [&wp](const typename H5::Complex& a) { return abssquare(a) / wp->N0; }); // hh
	std::transform(wp->g.begin(), wp->g.end(), gabs.begin(), [&wp,P0](const typename H5::Complex& a) { return abssquare(a) * P0 / wp->N; }); // gg

	for (size_t i = 0; i < D; ++i)
		this->setUB(i, P[i]);
	this->setLB(0);

	for (size_t i = 0; i < D; ++i)
	{
		double ginv = 1.0/gabs[mx_q[i]];

		this->w[i] = 1;
		this->alpha[i] = habs[i];
		this->sigma[i] = 1 + ginv;

		for (size_t j = 0; j < D; ++j)
			this->beta[i][j] = habs[j] * ginv;
		this->beta[i][mx_l[i]] += habs[mx_l[i]];
	}

	this->setPrecision(1e-3);
	this->useRelTol = false;
	this->disableReduction = true;
}


template <size_t D>
double
TIN_Dinkelbach<D>::obj(const vtype& x, const vtype& y) const
{
	double ret = 1;

	for (size_t i = 0; i < D; ++i) {
		double tmp = alpha[i] * x[i] + std::inner_product(x.begin(), x.end(), beta[i].begin(), sigma[i]);
		tmp /= std::inner_product(y.begin(), y.end(), beta[i].begin(), sigma[i]);
		if (w[i] != 1) {
			tmp = std::pow(tmp, w[i]);
		}
		ret *= tmp;
	}

	ret = std::log2(ret);
	ret -= mu * lambda * std::accumulate(y.begin(), y.end(), 0.0);

	return ret;
}


template <size_t D>
typename TIN_Dinkelbach<D>::H5::Result
TIN_Dinkelbach<D>::solve(const size_t)
{
	unsigned long long outerIter = 0;
	unsigned long long innerIter = 0;
	double cbv = 0;
	std::chrono::time_point<clk> tic = clk::now();

	lambda = 0;

	do
	{
		++outerIter;
		this->optimize();
		cbv = this->optval - lambda * Pc;

		double numer = 1;
		for (size_t i = 0; i < D; ++i)
		{
			double tmp = alpha[i] * this->xopt[i];
			tmp /= std::inner_product(this->xopt.begin(), this->xopt.end(), beta[i].begin(), sigma[i]);
			tmp += 1;
			if (w[i] != 1) {
				tmp = std::pow(tmp, w[i]);
			}
			numer *= tmp;
		}
		numer = std::log2(numer);

		double denom = mu * std::accumulate(this->xopt.begin(), this->xopt.end(), 0.0) + Pc;

		lambda = numer / denom;

		innerIter += this->iter;
	} while (std::abs(cbv) > outerTol);

	// output
	std::chrono::time_point<clk> toc = clk::now();

	typename H5::Result res;
	res.objective = lambda;
	res.iter = outerIter;
	res.lastUpdate = innerIter;
	res.runtime = std::chrono::duration<double>(toc - tic).count();

	std::copy_n(this->xopt.begin(), res.xopt_nc.size(), res.xopt_nc.begin());

	switch (this->status) {
	case BRB<D>::Status::Optimal:
		res.status = H5::SolverStatus::Optimal;
		break;

	case BRB<D>::Status::Infeasible:
		res.status = H5::SolverStatus::Infeasible;
		break;

	case BRB<D>::Status::Unsolved:
		res.status = H5::SolverStatus::Unsolved;
		break;
	}

	std::cout << std::endl << std::endl;
	std::cout << "lambda: " << lambda << std::endl;
	std::cout << "iter: " << outerIter << std::endl;
	std::cout << "total inner iter: " << innerIter << std::endl;
	std::cout << "Time: " << res.runtime << std::endl;

	return res;
}


#if 0
template <size_t D>
class Dinkelbach
{
	public:
		//using H5 = HDF<SNDDIM,SNDDIM>;
		using S = TIN_DM<D>;


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


	private:
		H5::WPPtr wp;

};

template <size_t D>
void
Dinkelbach<D>::solve(const size_t)
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

			cbv = s->optval - lambda * Pc;
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
#endif

} // end namespace mmp

using namespace mmp;
	

int
main(int argc, char *argv[])
{
	Sim<TIN_Dinkelbach<3>> s(argc, argv);
	s.run();
}
