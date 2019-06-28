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

#ifndef _HKSOLVER_H
#define _HKSOLVER_H

#include <memory>
#include <algorithm>
#include <iostream>
#include <sstream>
using std::cout; using std::endl;

#include "HK.h"
#include "HDF.h"

template <class HK = ::HK>
class HKSolver
{
public:
	HK s;

	using H5 = typename HK::H5;

	template <typename... Ts>
	HKSolver(typename H5::WPPtr wp, Ts&&... args) try
		: s(std::forward<Ts>(args)...)
	{
		s.init(wp);

		typename HK::vtypeC linC = {};
		std::fill_n(linC.begin(), HKDIM, 1.0);
		
		s.setObjective(linC, {});
	}
	catch (GRBException &e)
	{
		std::stringstream ss;
		ss << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")";
		throw std::runtime_error(ss.str());
	}

	virtual ~HKSolver() {};

	virtual std::string getName() const
	{
		return "hk";
	}

	size_t getDim()
	{
		return 1;
	}

	const HK* solve(const size_t)
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
