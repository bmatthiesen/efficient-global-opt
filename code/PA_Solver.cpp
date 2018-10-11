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

#include "PA_Solver.h"

#include <memory>

namespace pa {

userdata::userdata(const double* w)
	: model(GurobiEnv::getInstance())
{
	model.set(GRB_IntParam_LogToConsole, 0);

	R = model.addVars(NULL, NULL, w, NULL, NULL, Dim);
	model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

	for (size_t i = 0; i < Dim; i++)
		model.addConstr(R[i], GRB_LESS_EQUAL, 0.0);

	for (size_t i = 0; i < Dim; i++)
	{
		GRBLinExpr expr(R[i]);
		expr += GRBLinExpr(R[mx_l[i]]);

		model.addConstr(expr, GRB_LESS_EQUAL, 0.0);
	}

	model.update();

	constr = model.getConstrs();
}

userdata::~userdata()
{
	delete[] R;
	delete[] constr;
}

const unsigned short userdata::mx_l[] = {2,0,1};
const unsigned short userdata::mx_q[] = {1,2,0};



SNDSolver::SNDSolver(H5::WPPtr wp) try
	: s(), ud(w)
{
	auto P0 = std::pow(10, wp->P0/10.0);
	std::transform(wp->P.begin(), wp->P.end(), ud.P.begin(), [](const double p) { return std::pow(10, p/10.0);});

	std::transform(wp->h.begin(), wp->h.end(), ud.h.begin(), [&wp](const H5::Complex& a) { return abssquare(a) / wp->N0; });
	std::transform(wp->g.begin(), wp->g.end(), ud.g.begin(), [&wp,P0](const H5::Complex& a) { return abssquare(a) * P0 / wp->N; });

	{
		std::array<double,Dim> tmp;
		obj_g(ud.P, &tmp[0], &ud);
		ud.gsigma_max = std::accumulate(tmp.begin(), tmp.end(), 0.0);

		obj_g({0,0,0}, &tmp[0], &ud);
		ud.tmax = ud.gsigma_max - std::accumulate(tmp.begin(), tmp.end(), 0.0);
	}

	s.setUB({ud.P[0], ud.P[1], ud.P[2], ud.tmax});
	s.setObjective(obj);
	s.setInH(inH);
	s.setInG(inG);
	s.setUserData(&ud);
	s.setShift();
	s.setPrecision(1e-2);
}
catch (GRBException &e)
{
	std::stringstream ss;
	ss << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")";
	throw std::runtime_error(ss.str());
}

std::unique_ptr<SNDSolver::H5::Result>
SNDSolver::solve(const size_t)
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

	std::unique_ptr<H5::Result> res = std::make_unique<H5::Result>();

	res->objective = s.optval;
	std::copy_n(s.xopt.begin(), Dim, res->xopt_nc.begin());
	res->runtime = s.runtime;
	res->iter = s.iter;
	res->lastUpdate = s.lastUpdate;

	for (size_t i = 0; i < Dim; ++i)
		res->xopt_c[i] = ud.R[i].get(GRB_DoubleAttr_X);

	switch (s.status)
	{
		case PAstatus::OPTIMAL:
				res->status = H5::SolverStatus::Optimal;
				break;

		case PAstatus::MAXITER:
				res->status = H5::SolverStatus::Maxiter;
				break;

		case PAstatus::UNSOLVED:
				res->status = H5::SolverStatus::Unsolved;
				break;
	}

	return res;
}
	
	
void
objh(const vtype& x, double* res, void *user_data)
{
	userdata* ud = (userdata *) user_data;

	double dot = std::inner_product(ud->h.begin(), ud->h.end(), x.begin(), 1.0);
	
	for (size_t i = 0; i < Dim; i++)
		res[i] = 1 + dot / ud->g[ud->mx_q[i]];
}

void
obj_g(const vtype& x, double* res, void* user_data)
{
	double oh[Dim];
	objh(x, oh, user_data);

	obj_g(x, res, user_data, oh);
}

void
obj_g(const vtype&, double* res, void*, double* oh)
{
	for (size_t i = 0; i < Dim; i++)
		//res[i] = res[i+Dim] = log(oh[i]);
		res[i] = log(oh[i]);
}

basetype
obj(const vtype& x, void* user_data)
{
	userdata* ud = (userdata *) user_data;

	// compute h
	double oh[Dim];
	objh(x, oh, ud);

	// compute g
	double g[Dim];
	obj_g(x, g, ud, oh);

	// data for f
	std::array<double, Dim> snr;
	std::transform(ud->h.begin(), ud->h.end(), x.begin(), snr.begin(), [] (const double &a, const double &b) { return a*b; });

	// update constraints
	for (size_t i = 0; i < Dim; i++)
	{
		double gpsum = g[ud->mx_q[i]] + g[ud->mx_l[i]] + x[Dim] - ud->gsigma_max;

		double arg = oh[i] + snr[i];
		double foo = log(arg) + gpsum;
		ud->constr[i].set(GRB_DoubleAttr_RHS, foo);

		arg += snr[ud->mx_l[i]];
		foo = log(arg) + gpsum;
		ud->constr[i+Dim].set(GRB_DoubleAttr_RHS, foo);
	}

	ud->model.optimize();

	if (ud->model.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
		return -INF;
	else
		return ud->model.get(GRB_DoubleAttr_ObjVal);
}

basetype
g(const vtype& x, void* user_data)
{
	userdata* ud = (userdata *) user_data;
	const basetype* t = &x[Dim];

	std::array<double, Dim> g;
	obj_g(x, &g[0], ud);

	basetype res = std::accumulate(g.begin(), g.end(), *t - ud->gsigma_max);

	basetype tmp = *t - ud->tmax;
	if (tmp > res) res = tmp;

	for (size_t i = 0; i < Dim; i++)
	{
		tmp = x[i] - ud->P[i];
		if (tmp > res) res = tmp;
	}

	return res;
}

bool
inG(const vtype& x, void* user_data)
{
	userdata* ud = (userdata *) user_data;
	const basetype* t = &x[Dim];

	std::array<double, Dim> g;
	obj_g(x, &g[0], ud);

	return std::accumulate(g.begin(), g.end(), *t) <= ud->gsigma_max;
}

bool
inH(const vtype& x, void*)
{
	return x[0] >= 0 && x[1] >= 0 && x[2] >= 0 && x[3] >= 0;
}

}
