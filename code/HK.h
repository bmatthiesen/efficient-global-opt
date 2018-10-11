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

#ifndef _HK_H
#define _HK_H

#include "RR.h"
#include "HDF.h"

#define HKDIM 3

class HK : public RR<HKDIM+1,2*HKDIM,5*HKDIM+2>
{
	public:
		static const size_t dim = HKDIM;

		using _RR = RR<HKDIM+1,2*HKDIM,5*HKDIM+2>;
		using H5 = HDF<HKDIM+1,2*HKDIM,HKDIM>; // TODO

		static const unsigned short mx_l[];
		static const unsigned short mx_q[];

		/* ctor */
		HK(const BndSolvers bs = BndSolvers::Mosek);

		/* always call init */
		void init(H5::WPPtr wp);

		/* overrides */
		void optimize() override;

	private:
		bool initialized = false;

		template <typename T>
		static double abssquare(const T& a)
			{ return a.re * a.re + a.im * a.im; }
};

#endif

