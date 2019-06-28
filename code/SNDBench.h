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

#ifndef _SNDBENCH_H_
#define _SNDBENCH_H_

#include <set>
#include <iostream>
#include <algorithm>

#include "RR.h"
#include "HDFBench.h"

#define MAXUE 50


template <size_t numUE, size_t numNC, bool simple>
struct SNDBenchNumConstraints
{
	static constexpr size_t numC = numUE - numNC;
	static constexpr size_t cnt = (numC << numC) + (numNC << (numC+1)) - numUE;
};

template <size_t numUE, size_t numNC>
struct SNDBenchNumConstraints<numUE,numNC,true>
{
	static constexpr size_t numC = numUE - numNC;
	static constexpr size_t cnt = 3*numUE;
	static_assert(numC > 0, "This does not make sense (or is trivial)");
};

// NCDim = numNC
// CDim = numUE rates
// NumConstraints = numC * (2^numC - 1) + numNC * (2^(numC + 1) - 1)
template <size_t numUE, size_t numNC, bool simple = true>
class SNDBench : public RR<numNC, numUE, SNDBenchNumConstraints<numUE, numNC, simple>::cnt> // RR<NCDim, CDim, NumConstraints>
{
	static constexpr size_t numConstr = SNDBenchNumConstraints<numUE, numNC, simple>::cnt;

	static double N;
	static double P;

	using _RR = RR<numNC, numUE, numConstr>;

	public:
		using H5 = HDF<numNC, numUE, MAXUE>;
		using typename _RR::vtype;
		using typename _RR::vtypeC;

		static constexpr size_t oR = 0; // offset R

		// ctor
		SNDBench(const BndSolvers bs = BndSolvers::GurobiSep);

		// always call init
		void init(typename H5::WPPtr wp);

		/* overrides */
		void optimize() override;

	private:
		bool initialized = false;

		template <typename T>
		static double abssquare(const T& a)
			{ return a.re * a.re + a.im * a.im; }
};

template <class S>
auto powerset(const S& s);

template <typename T>
std::ostream& operator<<(std::ostream &out, const std::set<T>& s);

#include "bits/SNDBench.cpp"

#endif
