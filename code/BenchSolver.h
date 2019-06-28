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

#ifndef _BENCHSOLVER_H
#define _BENCHSOLVER_H

#include <memory>
#include <algorithm>
#include <iostream>
#include <sstream>
using std::cout; using std::endl;

#include "SNDBench.h"
#include "HDFBench.h"

template <size_t numUE, size_t numNC>
class BenchSolver
{
public:
	using _S = SNDBench<numUE, numNC>;
	using H5 = typename _S::H5;

	_S s;

	template <typename... Ts>
	BenchSolver(typename H5::WPPtr wp, Ts&&... args) try
		: s(std::forward<Ts>(args)...)
	{
		s.init(wp);

		typename _S::vtypeC linC = {};
		for (size_t i = 0; i < numUE; ++i)
			linC[s.oR + i] = 1.0;

		s.setObjective(linC, {});
	}
	catch (GRBException &e)
	{
		std::stringstream ss;
		ss << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")";
		throw std::runtime_error(ss.str());
	}

	virtual ~BenchSolver() {};

	virtual std::string getName() const
	{
		std::stringstream ss;
		ss << "sndBench-numUE_" << numUE << "-numNC_" << numNC;
		return ss.str();
	}

	size_t getDim()
	{
		return 1;
	}

	const _S* solve(const size_t)
	{
		try
		{
			s.optimize();
		}
		catch (GRBException &e)
		{
			std::stringstream ss;
			ss << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")";
			throw std::runtime_error(ss.str());
		}

		return &s;
	}

	void enableBackup(const std::string& filename, const std::chrono::seconds& interval)
	{
		s.enableBackup(filename, interval);
	}

	void removeBackup() const
	{
		s.removeBackup();
	}
};

#endif

