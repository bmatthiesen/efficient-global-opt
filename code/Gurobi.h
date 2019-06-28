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

#ifndef _GUROBI_H
#define _GUROBI_H

#include "gurobi_c++.h"

class GurobiEnv
{
	public:
		static size_t maxTries;
		static GRBEnv& getInstance();

		GurobiEnv() = delete;
		GurobiEnv(GurobiEnv const&) = delete;
		void operator=(GurobiEnv const&) = delete;

	private:
		static GRBEnv& _getInstance();
};

#endif
