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

#ifndef _PA_SOLVER_H
#define _PA_SOLVER_H

#include <sstream>
#include "gurobi_c++.h"

#include "Gurobi.h"

#include "PA.h"
#include "HDF.h"

namespace pa {

const size_t Dim = 3;

using vtype = PA<Dim+1>::vtype;

struct userdata
{
	static const unsigned short mx_l[Dim];
	static const unsigned short mx_q[Dim];

	userdata(const double* w);
	~userdata();

	// gurobi modell (auto init)
	GRBModel model;
	GRBVar *R;
	GRBConstr *constr;

	// user set variables
	std::array<double, Dim> h, g;
	double gsigma_max;
	vtype P;
	double tmax;
};

class SNDSolver
{
	public:
		using H5 = HDF<Dim,Dim>;

		SNDSolver(H5::WPPtr wp);

		std::string getName() const
		{
			return "afsndPA";
		}

		constexpr size_t getDim()
		{
			return 1;
		}

		std::unique_ptr<H5::Result> solve(const size_t);

		void enableBackup(const std::string&, const std::chrono::seconds&)
		{
			std::cerr << "enableBackup() not implemented" << std::endl;
		}

		void removeBackup() const
		{
			std::cerr << "removeBackup() not implemented" << std::endl;
		}
	private:
		PA<Dim+1> s;
		const double w[Dim] = {1, 1, 1};
		userdata ud;
};

const double INF = std::numeric_limits<double>::infinity();

inline double
log(double x)
{
	return x <= 0 ? -INF : std::log(x);
}

template <typename T>
constexpr double
abssquare(const T& a)
{
	return a.re * a.re + a.im * a.im;
}

void objh(const vtype& x, double* res, void *user_data);
void obj_g(const vtype& x, double* res, void* user_data);
void obj_g(const vtype& x, double* res, void*, double* oh);
basetype obj(const vtype& x, void* user_data);
basetype g(const vtype& x, void* user_data);
bool inG(const vtype& x, void* user_data);
bool inH(const vtype& x, void*);

}
#endif
