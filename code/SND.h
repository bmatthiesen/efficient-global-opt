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

#ifndef _SND_H
#define _SND_H

#include <bitset>

#include "RR.h"
#include "HDF.h"

#define SNDDIM 3

template <size_t c>
inline constexpr size_t
_bits()
{
	static_assert(c<=7, "out of range");
	return ((c>>2)&1) + ((c>>1)&1) + (c&1);
}
template <size_t c>
inline constexpr size_t bits = _bits<c>();


template <unsigned long long sndtype = 7>
class SND : public RR<SNDDIM,SNDDIM,SNDDIM+bits<sndtype>>
{
	public:
		static const size_t dim = SNDDIM;
	private:
		using _RR = RR<SNDDIM,SNDDIM,SNDDIM+bits<sndtype>>;
		using H5 = HDF<SNDDIM,SNDDIM>;
	public:
		using typename _RR::vtype;
		using typename _RR::vtypeC;

		static const unsigned short mx_l[];
		static const unsigned short mx_q[];

		/* ctors */
		SND(const BndSolvers bs = BndSolvers::GurobiSep);

		/* always call init */
		void init(H5::WPPtr wp, const bool useTighterBound = true);

		/* overrides */
		void optimize() override;
			
	private:
		bool initialized = false;
		std::bitset<dim> useSND;

		template <typename T>
		static double abssquare(const T& a)
			{ return a.re * a.re + a.im * a.im; }

};

#include "bits/SND.cpp"

#endif
