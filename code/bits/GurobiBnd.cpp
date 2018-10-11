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

// included from RR.h

template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::GurobiBnd::GurobiBnd(const size_t maxLogTerms)
	: 	Grb(GurobiEnv::getInstance()),
		varsC(Grb.addVars(CDim))
{
	logterms.reserve(maxLogTerms);

	Grb.set(GRB_IntParam_LogToConsole, 0);
	Grb.setObjective(GRBLinExpr(0), GRB_MINIMIZE);
	
	varT = Grb.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);

	Grb.update();
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::GurobiBnd::putvarbound_NC(const size_t varIdx, const basetype lb, const basetype ub)
{
	if (NCubSgnU[varIdx] == 0)
		throw std::runtime_error("foo");
	else if (NCubSgnU[varIdx] > 0) // == 0 is probably an error TODO check and remove previous test
		NCub[varIdx] = ub;
	else
		NCub[varIdx] = lb;
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::GurobiBnd::putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val)
{
	if (val != 0)
	{
		if (NCubSgnU[varIdx] == 0)
			NCubSgnU[varIdx] = val < 0 ? 1 : -1;
		else if (NCubSgnU[varIdx] > 0 && val > 0.0)
		{
			std::stringstream ss;
			ss << "linNC[" << varIdx << "] must must be <= 0.0";
			throw std::runtime_error(ERR(ss.str()));
		}
		else if (NCubSgnU[varIdx] < 0 && val < 0.0)
		{
			std::stringstream ss;
			ss << "linNC[" << varIdx << "] must must be >= 0.0";
			throw std::runtime_error(ERR(ss.str()));
		}
	}

	constr_linNC[constrIdx][varIdx] = val;
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::GurobiBnd::setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd)
{
	for (size_t i = 0; i < NCDim; ++i)
	{
		if (NCubSgnU[i] == 0)
		{
			if (linNC[i] != 0)
				NCubSgnU[i] = linNC[i] < 0 ? 1 : -1;
		}
		else if (linNC[i] * NCubSgnU[i] > 0)
		{
			std::stringstream ss;
			ss << "linNC must have signs " << NCubSgnU;
			throw std::runtime_error(ERR(ss.str()));
		}
	}

	constr_linNC[constrIdx] = linNC;
	constr_consts[constrIdx] = c;

	GRBLinExpr lhs;
	lhs.addTerms(linC.data(), varsC.get(), linC.size());
	
	char sense;
	switch (bnd)
	{
		case Boundtype::Upper:
			sense = GRB_LESS_EQUAL;
			break;

		default:
			throw std::runtime_error("boundtype not supported");
			break;
	}

	constr_linC[constrIdx] = std::move(Grb.addConstr(lhs, sense, c));
}



template <size_t NCDim, size_t CDim, size_t NumConstraints>
size_t
RR<NCDim,CDim,NumConstraints>::GurobiBnd::setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign)
{
	short factor = negativeSign ? -1 : 1;

	for (size_t i = 0; i < NCDim; ++i)
	{
		if (NCubSgnU[i] == 0)
		{
			if (linNC[i] != 0)
				NCubSgnU[i] = linNC[i] > 0 ? factor : -1*factor;
		}
		else if (linNC[i] * NCubSgnU[i] * factor < 0)
		{
			std::stringstream ss;
			ss << "linNC must have signs " << NCubSgnU;
			throw std::runtime_error(ERR(ss.str()));
		}
	}

	if (std::any_of(linC.begin(), linC.end(), [] (basetype e) { return e != 0.0; }))
		throw std::runtime_error(ERR("linC must be == 0.0"));

	logterms.emplace_back(linNC, c, constrIdx, negativeSign);

	return logterms.size()-1;
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::GurobiBnd::updateModel()
{
	auto constr = constr_consts;

	for (size_t i = 0; i < constr_linNC.size(); ++i)
		constr[i] -= std::inner_product(constr_linNC[i].begin(), constr_linNC[i].end(), NCub.begin(), 0.0);

	for (auto& l : logterms)
	{
		auto val = std::log(std::inner_product(l.linNC.begin(), l.linNC.end(), NCub.begin(), l.c));

		if (l.negativeSign)
			constr[l.constrIdx] -= val;
		else
			constr[l.constrIdx] += val;
	}


	for (size_t i = 0; i < constr_linC.size(); ++i)
		constr_linC[i].set(GRB_DoubleAttr_RHS, constr[i]);

	Grb.update();
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
bool
RR<NCDim,CDim,NumConstraints>::GurobiBnd::solve()
{
	updateModel();
	Grb.optimize();

	return Grb.get(GRB_IntAttr_Status) == GRB_OPTIMAL;
}
