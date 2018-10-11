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

// include from RR.h

template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::GurobiObj::GurobiObj()
	: 	Grb(GurobiEnv::getInstance()),
		varsC(Grb.addVars(CDim))
{
	Grb.set(GRB_IntParam_LogToConsole, 0);
	Grb.setObjective(GRBLinExpr(0), GRB_MAXIMIZE);
	
	Grb.update();
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::GurobiObj::setConstr(const size_t constrIdx, const vtypeC& linC, const basetype c, const Boundtype bnd)
{
	char sense;
	switch (bnd)
	{
		case Boundtype::Upper:
			sense = GRB_LESS_EQUAL;
			break;

		case Boundtype::Lower:
			sense = GRB_GREATER_EQUAL;
			break;

		case Boundtype::Fixed:
			sense = GRB_EQUAL;
			break;

		default:
		case Boundtype::Undefined:
			throw std::runtime_error("boundtype not supported");
			break;
	}

	GRBLinExpr lhs;
	lhs.addTerms(linC.data(), varsC.get(), linC.size());

	constr_linC[constrIdx] = std::move(Grb.addConstr(lhs, sense, c));
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
size_t
RR<NCDim,CDim,NumConstraints>::GurobiObj::setLogConstr(const size_t, const vtypeC&, const basetype, const bool)
{
	throw std::runtime_error(ERR("not supported by linear solver"));
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::GurobiObj::updateLogConstant(const size_t, const basetype)
{
	throw std::runtime_error(ERR("not supported by linear solver"));
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::GurobiObj::putobj(const vtypeC& numerLinC, const vtypeC& denomLinC)
{
	if (!std::all_of(denomLinC.begin(), denomLinC.end(), [] (basetype e) { return e == 0.0; }))
		throw std::runtime_error(ERR("denomLinC is only partially implemented. Please consult source code comments near this message for more information."));
	else
	{
		for (size_t i = 0; i < CDim; ++i)
			varsC[i].set(GRB_DoubleAttr_Obj, numerLinC[i]);
	}
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
inline void
RR<NCDim,CDim,NumConstraints>::GurobiObj::putobjConst(const basetype numer, const basetype denom)
{
	numerConst = numer;
	denomConst = denom;
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
bool
RR<NCDim,CDim,NumConstraints>::GurobiObj::solve()
{
	updateModel();
	Grb.optimize();

	return Grb.get(GRB_IntAttr_Status) == GRB_OPTIMAL;
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::GurobiObj::getopt_C(vtypeC& x) const
{
	for (size_t i = 0; i < x.size(); ++i)
		x[i] = varsC[i].get(GRB_DoubleAttr_X);
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
basetype
RR<NCDim,CDim,NumConstraints>::GurobiObj::getobj() const
{
	basetype sol = Grb.get(GRB_DoubleAttr_ObjVal);

	return (sol + numerConst) / denomConst;
}
