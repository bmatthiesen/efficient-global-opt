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

template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::MosekSepBnd::MosekSepBnd(const size_t maxLogTerms)
	: MosekBase(CDim+1, CDim, 0, NumConstraints+1, 0, "bnd"),
		vOffsetT(CDim),
		constr_consts({}),
		constr_linNC({}),
		NCub({})
{
	logterms.reserve(maxLogTerms);

	mskcall(MSK_putvarbound(task, vOffsetT, MSK_BK_FR, -MSK_INFINITY, +MSK_INFINITY));
	MosekBase::putobj(vOffsetT, 1.0);
	//mskcall(MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr));
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekSepBnd::setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd)
{
	if (std::any_of(linNC.begin(), linNC.end(), [] (basetype e) { return e > 0.0; }))
		throw std::runtime_error(ERR("linNC must be <= 0.0"));

	if (bnd != Boundtype::Upper)
		throw std::runtime_error(ERR("unsupported boundtype"));

	constr_linNC[constrIdx] = linNC;
	constr_consts[constrIdx] = c;

	MosekBase::setConstr(constrIdx, linC, {}, c, bnd);
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekSepBnd::putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val)
{
	if (val > 0.0)
		throw std::runtime_error(ERR("linNC must be <= 0.0"));

	constr_linNC[constrIdx][varIdx] = val;
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
size_t
RR<NCDim,CDim,NumConstraints>::MosekSepBnd::setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign)
{
	if (negativeSign)
		throw std::runtime_error(ERR("not implemented"));

	if (std::any_of(linNC.begin(), linNC.end(), [] (basetype e) { return e < 0.0; }))
		throw std::runtime_error(ERR("linNC must be >= 0.0"));

	if (std::any_of(linC.begin(), linC.end(), [] (basetype e) { return e != 0.0; }))
		throw std::runtime_error(ERR("linC must be == 0.0"));

	logterms.emplace_back(linNC, c, constrIdx);

	return logterms.size()-1;
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekSepBnd::putconbound(const size_t constrIdx, const basetype ub)
{
	constr_consts[constrIdx] = ub;
	MosekBase::putconbound(constrIdx, ub);
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekSepBnd::updateModel()
{
	auto constr = constr_consts;
	std::vector<bool> update(constr_consts.size(), false);

	for (auto& l : logterms)
	{
		auto val = std::log(std::inner_product(l.linNC.begin(), l.linNC.end(), NCub.begin(), l.c));

		if (val == 0.0)
			continue;

		constr[l.constrIdx] += val;
		update[l.constrIdx] = true;
	}

	for (size_t i = 0; i < update.size(); ++i)
	{
		if (update[i])
			MosekBase::putconbound(i, constr[i]);
	}

	MosekBase::updateModel();
}

