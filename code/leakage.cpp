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
#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include <limits>
#include "gurobi_c++.h"

#include "SIT.h"

using std::cout;
using std::cerr;
using std::endl;
using std::pow;

class Ex : public SIT<2>
{
	public:
		Ex(const double P, const double L, const double eps)
			: SIT(), leakage(L*pow(2, eps)), qos(6*pow(2, -eps)), env(), Gbnd(env)
		{
#ifdef REDUCTION
			disableReduction = true;
#endif
			gamma0 = -std::numeric_limits<double>::infinity();

			buildBnd();

			setPrecision(1e-6, 1e-6);
			setLB(0);
			setUB(P);
		}

		~Ex()
		{
			delete[] Gbnd_x;
		}

	private:
		const double leakage;
		const double qos;
		const double c[2] { -1, 0 }; // P, Q

#ifdef REDUCTION
		void reduction(const PBox&, RBox&) override final {};
#endif
		void bound(RBox& red) override final;

		basetype obj(const vtype& P, bool update) override final;
		bool isFeasible(const vtype& P) override final;

		// Gurobi
		GRBEnv env;
		GRBModel Gbnd;
		GRBVar *Gbnd_x;
		GRBConstr Gbnd_gamma;

		void buildBnd(void);
		void setGamma(basetype g) override final;
};


void
Ex::buildBnd(void)
{
	Gbnd.set(GRB_IntParam_LogToConsole, 0);
	
	Gbnd_x = Gbnd.addVars(2);
	GRBVar &p = Gbnd_x[0];
	GRBVar &q = Gbnd_x[1];

	Gbnd.setObjective(.25*p*p + .25*q*q + .5*p*q + .5*p + q, GRB_MINIMIZE);

	Gbnd.addConstr(p + q, GRB_GREATER_EQUAL, qos);

	Gbnd_gamma = Gbnd.addConstr(c[0]*p + c[1]*q, GRB_GREATER_EQUAL, 0, "gamma");

	Gbnd.update();
}


void
Ex::setGamma(basetype g)
{
	SIT::setGamma(g);

	Gbnd_gamma.set(GRB_DoubleAttr_RHS, g);
}


void
Ex::bound(RBox& red)
{
	for (size_t i = 0; i < 2; i++)
	{
		Gbnd_x[i].set(GRB_DoubleAttr_LB, red.lb(i));
		Gbnd_x[i].set(GRB_DoubleAttr_UB, red.ub(i));
	}

	Gbnd.optimize();

	if (Gbnd.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		for (size_t i = 0; i < 2; i++)
			red.xk(i) = Gbnd_x[i].get(GRB_DoubleAttr_X);

		red.yk() = &red.ub();

		red.beta = Gbnd.get(GRB_DoubleAttr_ObjVal);
		red.beta += 1 - leakage - .25 * (pow(red.yk(0),2) + pow(red.yk(1),2));
	}
	else
	{
		red.beta = std::numeric_limits<basetype>::infinity();
	}
}

basetype
Ex::obj(const vtype& P, bool)
{
	return c[0]*P[0] + c[1]*P[1];
}


bool
Ex::isFeasible(const vtype& P)
{
	const double g2 = .5*P[0]*P[1] + .5*P[0] + P[1] + 1 - leakage;
	return (g2 <= 0.0);
}

int main(int argc, char *argv[])
{
	if (argc < 2 || argc > 3)
	{
		cerr << argv[0] << ": Wrong number of arguments: " << argv[0] << " leakage [epsilon]" << endl;
		return 2;
	}

	double leakage = std::strtod(argv[1], nullptr);
	double eps = argc == 3 ? std::strtod(argv[2],nullptr) : 0;

	Ex sit(5, leakage, eps);

	sit.optimize();

	return 0;
}
