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
std::ostream& operator<< (std::ostream &out, const typename RR<NCDim,CDim,NumConstraints>::Boundtype &point)
{
	switch (point)
	{
		case RR<NCDim,CDim,NumConstraints>::Boundtype::Undefined:
			out << "Undefined";
			break;

		case RR<NCDim,CDim,NumConstraints>::Boundtype::Lower:
			out << "Lower";
			break;

		case RR<NCDim,CDim,NumConstraints>::Boundtype::Upper:
			out << "Upper";
			break;

		case RR<NCDim,CDim,NumConstraints>::Boundtype::Fixed:
			out << "Fixed";
			break;

		default:
			out << "ERROR";
			break;
	}

	return out;
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::RR(const size_t maxLogterms, const BndSolvers bndS, const ObjSolvers objS)
	:	SIT<NCDim>(),
		maxLog(maxLogterms)
{
#ifdef REDUCTION
	red = std::make_unique<MosekRed>(maxLogterms);
#endif

	gamma0 = 0;

	switch (bndS) {
	case BndSolvers::Mosek:
		bnd = std::make_unique<MosekBnd>(maxLogterms); // vars: NC + C + t
		break;
	case BndSolvers::MosekSep:
		bnd = std::make_unique<MosekSepBnd>(maxLogterms); // vars: NC + C + t
		break;
	case BndSolvers::GurobiSep:
		bnd = std::make_unique<GurobiBnd>(maxLogterms); // vars: NC + C + t
		break;
	}

	switch (objS) {
	case ObjSolvers::Mosek:
		objs = std::make_unique<MosekObj>(maxLogterms);
		break;
	case ObjSolvers::Gurobi:
		objs = std::make_unique<GurobiObj>();
		break;
	}
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::setVarbounds_NC(const size_t varIdx, const basetype lb, const basetype ub)
{
	/* error checks */
	if (varIdx < 0 || varIdx >= NCDim)
		throw std::out_of_range(ERR("varIdx"));

	if (std::isnan(lb))
		throw std::logic_error(ERR("lb not defined"));

	if (std::isnan(ub))
		throw std::logic_error(ERR("ub not defined"));

	/* SIT */
	this->setLB(varIdx, lb);
	this->setUB(varIdx, ub);
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::setVarbounds_C(const size_t varIdx, const basetype lb, const basetype ub)
{
	/* error checks */
	if (varIdx < 0 || varIdx >= CDim)
		throw std::out_of_range(ERR("varIdx"));

	/* error checks */
	if (std::isnan(lb))
		throw std::logic_error(ERR("lb not defined"));

	if (lb < 0)
		throw std::logic_error(ERR("lb < 0"));

	if (!std::isnan(ub))
	{
		if (lb >= ub)
			throw std::logic_error(ERR("lb >= ub"));

		/* put range bound */
		bnd->putvarbound_C(varIdx, lb, ub);
#ifdef REDUCTION
		red->putvarbound_C(varIdx, lb, ub);
#endif
		objs->putvarbound_C(varIdx, lb, ub);
	}
	else
	{
		/* put lower bound */
		bnd->putvarbound_C(varIdx, lb);
#ifdef REDUCTION
		red->putvarbound_C(varIdx, lb);
#endif
		objs->putvarbound_C(varIdx, lb);
	}

	setStatus(SIT<NCDim>::Status::Unsolved);
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::setObjective(const vtypeC& numerLinC, const vtype& numerLinNC, const basetype numerConst)
{
	obj_numerLinNC = numerLinNC;
	obj_numerLinC = numerLinC;
	obj_numerConst = numerConst;

	obj_denomLinNC = {};
	obj_denomLinC = {};
	obj_denomConst = 0.0;
	obj_denomPresent = false;

	setObjective();
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::setObjective(const vtypeC& numerLinC, const vtype& numerLinNC, const vtypeC& denomLinC, const vtype& denomLinNC, const basetype numerConst, const basetype denomConst)
{
	obj_numerLinNC = numerLinNC;
	obj_numerLinC = numerLinC;
	obj_numerConst = numerConst;

	obj_denomLinNC = denomLinNC;
	obj_denomLinC = denomLinC;
	obj_denomConst = denomConst;

	auto allZero = [] (const auto& v) { return std::all_of(v.begin(), v.end(), [] (basetype e) { return e == 0.0; }); };
	obj_denomPresent = !(obj_denomConst == 0.0 && allZero(obj_denomLinNC) && allZero(obj_denomLinC));

	setObjective();
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::setObjective()
{
	setStatus(SIT<NCDim>::Status::Unsolved);

	/* obj */
	objs->putobj(obj_numerLinC, obj_denomLinC);

	/* pass to solver */
	vtypeC linC;
	vtype linNC;
	if (obj_denomPresent)
	{
		std::transform(obj_numerLinNC.begin(), obj_numerLinNC.end(), obj_denomLinNC.begin(), linNC.begin(), [this] (const auto numer, const auto denom) { return gamma*denom - numer; });
		std::transform(obj_numerLinC.begin(), obj_numerLinC.end(), obj_denomLinC.begin(), linC.begin(), [this] (const auto numer, const auto denom) { return gamma*denom - numer; });
	}
	else
	{
		std::transform(obj_numerLinNC.begin(), obj_numerLinNC.end(), linNC.begin(), [this] (const auto numer) { return -numer; });
		std::transform(obj_numerLinC.begin(), obj_numerLinC.end(), linC.begin(), [this] (const auto numer) { return -numer; });
	}

#ifdef REDUCTION
	red->setConstr(cOffsetObj, linC, linNC, obj_numerConst - gamma * obj_denomConst, Boundtype::Upper);
#endif
	bnd->setConstr(cOffsetObj, linC, linNC, obj_numerConst - gamma * obj_denomConst, Boundtype::Upper);
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::setConstraint(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd)
{
	/* error checks */
	if (constrIdx < 0 || constrIdx >= NumConstraints)
		throw std::out_of_range(ERR("constrIdx"));

	setStatus(SIT<NCDim>::Status::Unsolved);

	/* save */
	constrConsts[constrIdx] = c;

	//if (std::any_of(linNC.begin(), linNC.end(), [](basetype e) { return e != 0.0; }))
		constr_linNC[constrIdx] = linNC;
	//else if (!constr_linNC[constrIdx].empty())
		//constr_linNC[constrIdx].clear();

	/* pass to solver */
	objs->setConstr(constrIdx, linC, c, bnd);
#ifdef REDUCTION
	red->setConstr(cOffsetC + constrIdx, linC, linNC, c, bnd);
#endif
	this->bnd->setConstr(cOffsetC + constrIdx, linC, linNC, c, bnd);
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::setLogterm(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const vtype& dcLinNC, const basetype c, const basetype dcC, const bool negativeSign)
{
	/* error checks */
	if (constrIdx < 0 || constrIdx >= NumConstraints)
		throw std::out_of_range(ERR("constrIdx"));

	for (auto e : dcLinNC)
		if (e < 0.0)
			throw std::logic_error(ERR("Negative coefficient in dcLinNC not implemented"));

	setStatus(SIT<NCDim>::Status::Unsolved);


	Logterm l(linC, linNC, dcLinNC, c, dcC, negativeSign);

	/* pass */
	if (!l.allzero_dcLinNC)
		// include t in bnd
		bnd->putaij_T(cOffsetC + constrIdx, -1.0);

	if (!l.allzero_linC)
		l.obj_logIdx = objs->setLogConstr(constrIdx, linC, c, negativeSign);

	if (!(l.allzero_linC && l.allzero_linNC))
	{
#ifdef REDUCTION
		red->setLogConstr(cOffsetC + constrIdx, linC, linNC, c, negativeSign);
#endif
		bnd->setLogConstr(cOffsetC + constrIdx, linC, linNC, c, negativeSign);
	}

	// die aux constraints erstellt solver
	/* store logterm */
	logterms[constrIdx].push_back(std::move(l));
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::bound(RBox& red)
{
	for (size_t i = 0; i < NCDim; ++i)
		bnd->putvarbound_NC(i, red.lb(i), red.ub(i));

	vtype* b = nullptr;

	for (size_t i = 0; i < NumConstraints; ++i)
	{
		basetype geval(constrConsts[i]);

		for (auto& l : logterms[i])
		{
			basetype tmp;

			if (!l.allzero_dcLinNC)
			{
				if (b == nullptr)
					b = l.negativeSign ? &red.ub() : &red.lb();

				tmp = std::log(std::inner_product(b->begin(), b->end(), l.dcLinNC.begin(), l.dcC));
			}
			else
				tmp = std::log(l.dcC);

			if (l.negativeSign)
				geval += tmp;
			else
				geval -= tmp;
		}

		bnd->putconbound(cOffsetC + i, geval);
	}

	if (bnd->solve())
	{
		red.beta = bnd->getobj();
		bnd->getopt_NC(red.xk());

		if (b == nullptr)
			throw std::runtime_error(ERR("cannot set red.yk"));

		red.yk() = b;
	}
	else
		red.beta = std::numeric_limits<basetype>::infinity();
}

#ifdef REDUCTION
template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::reduction(const PBox& P, RBox& _red)
{
	// initialize
	for (size_t i = 0; i < NCDim; ++i)
	{
		// set box for P
		red->putvarbound_NC(i, P.lb[i], P.ub[i]);

		// set objective to zero
		red->putobj_NC(i, 0.0);
	}

	// set bounds
	for (size_t i = 0; i < NumConstraints; ++i)
	{
		basetype geval(epsilon + constrConsts[i]);

		for (auto& l : logterms[i])
		{
			basetype tmp;

			if (!l.allzero_dcLinNC)
			{
				const vtype& b = l.negativeSign ? P.ub : P.lb;
				tmp = std::log(std::inner_product(b.begin(), b.end(), l.dcLinNC.begin(), l.dcC));
			}
			else
				tmp = std::log(l.dcC);

			if (l.negativeSign)
				geval += tmp;
			else
				geval -= tmp;
		}

		red->putconbound(cOffsetC + i, geval);
	}

	// reduction
	for (size_t i = 0; i < NCDim; ++i)
	{
		// update objective
		if (i > 0)
			red->putobj_NC(i-1, 0.0);

		red->putobj_NC(i, 1.0);

		// p'_i
		red->setMinimize();
		_red.lb(i) = red->solve() ? red->getobj() : P.lb[i];

		// q'_i
		red->setMaximize();
		_red.ub(i) = red->solve() ? red->getobj() : P.ub[i];
	}
}
#endif


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::updateObj(const vtype& P)
{
	// objective
	basetype numerConst = std::inner_product(P.begin(), P.end(), obj_numerLinNC.begin(), obj_numerConst);
	basetype denomConst = 1.0;
	if (obj_denomPresent)
		denomConst =  std::inner_product(P.begin(), P.end(), obj_denomLinNC.begin(), obj_denomConst);

	objs->putobjConst(numerConst, denomConst);

	// constraints
	for (size_t i = 0; i < NumConstraints; ++i)
	{
		basetype eval(constrConsts[i]);

		//if (!constr_linNC[i].empty())
			eval -= std::inner_product(P.begin(), P.end(), constr_linNC[i].begin(), 0.0);

		for (auto& l : logterms[i])
		{
			basetype feval = std::inner_product(P.begin(), P.end(), l.linNC.begin(), l.c);
			basetype geval = std::inner_product(P.begin(), P.end(), l.dcLinNC.begin(), l.dcC);

			basetype tmp;
			if (l.allzero_linC)
				tmp = std::log(feval / geval);
			else
			{
				objs->updateLogConstant(l.obj_logIdx, feval);
				tmp = -std::log(geval);
			}

			if (l.negativeSign)
				eval -= tmp;
			else
				eval += tmp;
		}

		objs->putconbound(i, eval);
	}
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
basetype
RR<NCDim,CDim,NumConstraints>::obj(const vtype& P, vtypeC* R, bool update)
{
	bool ret;

	if (update)
	{
		updateObj(P);
		ret = objs->solve();
	}
	else
		ret = true;

	if (ret)
	{
		basetype sol = objs->getobj();

		if (R != nullptr)
			objs->getopt_C(*R);

		return sol;
	}
	else
		return std::numeric_limits<basetype>::infinity();
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
bool
RR<NCDim,CDim,NumConstraints>::isFeasible(const vtype& P) // check g(P) <= 0.0 (feasibility problem)
{
	updateObj(P);
	return objs->solve();
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::setGamma(const basetype g)
{
	SIT<NCDim>::setGamma(g);

	if (obj_denomPresent)
	{
		for (size_t i = 0; i < NCDim; ++i)
		{
			basetype n(0.0), d(0.0);

			//if (!obj_numerLinNC.empty())
				n = obj_numerLinNC[i];

			//if (!obj_denomLinNC.empty())
				d = obj_denomLinNC[i];

			auto tmp = g * d - n;

#ifdef REDUCTION
			red->putaij_NC(cOffsetObj, i, tmp);
#endif
			bnd->putaij_NC(cOffsetObj, i, tmp);
		}

		for (size_t i = 0; i < CDim; ++i)
		{
			auto tmp = g * obj_denomLinC[i] - obj_numerLinC[i];

#ifdef REDUCTION
			red->putaij_C(cOffsetObj, i, tmp);
#endif
			bnd->putaij_C(cOffsetObj, i, tmp);
		}

#ifdef REDUCTION
		red->putconbound(cOffsetObj, obj_numerConst - g * obj_denomConst);
#endif
		bnd->putconbound(cOffsetObj, obj_numerConst - g * obj_denomConst);
	}
	else
	{
#ifdef REDUCTION
		red->putconbound(cOffsetObj, obj_numerConst - g);
#endif
		bnd->putconbound(cOffsetObj, obj_numerConst - g);
	}
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::checkModel()
{
	// check sign of logterms
	short ykSign = 0;

	for (auto& v : logterms)
	{
		for (auto& l : v)
		{
			if (!l.allzero_dcLinNC)
			{
				if (l.negativeSign)
				{
					if (ykSign == 0)
						ykSign = -1;
					else if (ykSign > 0)
						throw std::runtime_error(ERR("logterm sign wrong"));
				}
				else
				{
					if (ykSign == 0)
						ykSign = 1;
					else if (ykSign < 0)
						throw std::runtime_error(ERR("logterm sign wrong"));
				}
			}
		}
	}
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::optimize()
{
	checkModel();
#ifdef REDUCTION
	red->updateModel();
#endif
	bnd->updateModel();
	objs->updateModel();

	disablePrintResult = true;
	SIT<NCDim>::optimize();
	disablePrintResult = false;

	optval = obj(xopt, &_optC); // set optC

	if (output)
		printResult(); // optval in nats
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::printResult() const
{
	if (disablePrintResult)
		return;

	std::cout << "Status: " << statusStr << std::endl;
	std::cout << "Optval: " << optval << std::endl;

	std::cout << "x*_C: [";
	std::for_each(_optC.begin(), _optC.end(), [] (const basetype &a) { std::cout << " " << a; });
	std::cout << " ]" << std::endl;

	std::cout << "x*_NC: [";
	std::for_each(xopt.begin(), xopt.end(), [] (const basetype &a) { std::cout << " " << a; });
	std::cout << " ]" << std::endl;

	std::cout << "Precision: eta = " << eta << ", epsilon = " << epsilon << std::endl;

	std::cout << "Iter: " << iter << std::endl;
	std::cout << "Solution found in iter: " << lastUpdate << std::endl;

	std::cout << "Runtime: " << runtime << " sec" << std::endl;
}
