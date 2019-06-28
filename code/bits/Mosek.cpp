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

#if DEBUG
	#include <execinfo.h>
#endif

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::putvarbound(const size_t varIdx, const basetype lb, const basetype ub)
{
#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << ": putvarbound:" << lb << " < " << varIdx << " < " << ub << " : " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

	mskcall(MSK_putvarbound(task, varIdx, MSK_BK_RA, lb, ub));
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::putvarbound(const size_t varIdx, const basetype lb)
{
#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << ": putvarbound:" << lb << " < " << varIdx << " < INF : " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

	mskcall(MSK_putvarbound(task, varIdx, MSK_BK_LO, lb, +MSK_INFINITY));
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::putconbound(const size_t constrIdx, const basetype val)
{
#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << ": putconbound(" << constrIdx << ", " << val << ") / " << btypes[constrIdx] << ": " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

	switch (btypes[constrIdx])
	{
		case Boundtype::Upper:
			mskcall(MSK_putconbound(task, constrIdx, MSK_BK_UP, -MSK_INFINITY, val));
			break;

		case Boundtype::Lower:
			mskcall(MSK_putconbound(task, constrIdx, MSK_BK_LO, val, MSK_INFINITY));
			break;

		case Boundtype::Fixed:
			mskcall(MSK_putconbound(task, constrIdx, MSK_BK_FX, val, val));
			break;

		case Boundtype::Undefined:
		default:
			throw std::runtime_error("boundtype not set");
			break;
	}
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::putobj(const size_t varIdx, const basetype val)
{
#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << ": putcj(" << varIdx << ", " << val << ") : " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

	mskcall(MSK_putcj(task, varIdx, val));
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
basetype
RR<NCDim,CDim,NumConstraints>::MosekBase::getobj() const
{
	MSKrealt obj;
	mskcall(MSK_getprimalobj(task, MSK_SOL_ITR, &obj));

	return obj;
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::getopt(basetype* x, const size_t offset, const size_t len) const
{
	mskcall(MSK_getxxslice(task, MSK_SOL_ITR, offset, offset+len, x));
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
bool
RR<NCDim,CDim,NumConstraints>::MosekBase::solve()
{
	updateModel();

	try {
		mskcall(MSK_optimize(task));
	}
	catch (MSK_error& e) {
		if (e.rescode < 10000)
			throw;
#if DEBUG>3
		else
			std::cerr << e.what() << std::endl;
#endif
	}

	MSKsolstae solsta;
	mskcall(MSK_getsolsta(task, MSK_SOL_ITR, &solsta));

#if DEBUG > 4
	char str[2048];
	MSK_solstatostr(task, solsta, str);
	std::cout << name << ": sol status = " << str << std::endl;
#endif

	return solsta == MSK_SOL_STA_OPTIMAL;
}

#ifdef REDUCTION
template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekRed::setMinimize()
{
#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << ": Minimize : " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

	mskcall(MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE));
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekRed::setMaximize()
{
#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << ": Maximize : " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

	mskcall(MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE));
}
#endif

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::putaij(const size_t constrIdx, const size_t varIdx, const basetype val)
{
#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << ": putaij(" << constrIdx << ", " << varIdx << ", " << val << "): " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

	mskcall(MSK_putaij(task, constrIdx, varIdx, val));
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd)
{
	assert(linC.size() == dimC);
	assert(dimNC == 0 || linNC.size() == dimNC);

	for (size_t j = 0; j < dimC; j++)
	{
		if (linC[j] != 0.0)
			putaij(constrIdx, vOffsetC + j, linC[j]);
	}

	for (size_t j = 0; j < dimNC; j++)
	{
		if (linNC[j] != 0.0)
			putaij(constrIdx, vOffsetNC + j, linNC[j]);
	}

	assert(constrIdx < btypes.size());
	btypes[constrIdx] = bnd;
	putconbound(constrIdx, c);
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::MosekBase::MosekBase(const size_t numVars, const size_t Cdim, const size_t NCdim, const size_t numConstraints, const size_t maxLogTerms, const char* name) try
	:
		task(nullptr), btypes(numConstraints, Boundtype::Undefined),
		numVars(numVars), dimNC(NCdim), dimC(Cdim), logIdx(-1), numConstr(numConstraints), maxLog(maxLogTerms),
		vOffsetNC(0), vOffsetC(NCdim), 
		MSKsch(nullptr), update_MSKsch(false),
		oprc(new int[maxLogTerms]), opric(new int[maxLogTerms]), oprjc(new int[maxLogTerms]),
		oprfc(new double[maxLogTerms]), oprgc(new double[maxLogTerms]), oprhc(new double[maxLogTerms]),
		name(name)
{
	mskcall(MSK_maketask(MosekEnv::getEnv(), numConstraints+maxLogTerms, numVars+maxLogTerms, &task));

#if DEBUG > 5
	mskcall(MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr));
#endif

	mskcall(MSK_appendvars(task, numVars));
	mskcall(MSK_appendcons(task, numConstraints));
	mskcall(MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE));
	mskcall(MSK_putintparam(task, MSK_IPAR_INTPNT_MULTI_THREAD, MSK_OFF));
	mskcall(MSK_puttaskname(task, name));

}
catch (...)
{
	destroy();
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::MosekBase::~MosekBase() noexcept
{
	destroy();
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::destroy() noexcept
{
	if (MSKsch != nullptr)
		MSK_scend(task, &MSKsch);

	if (task != nullptr)
		MSK_deletetask(&task);
}


#ifdef REDUCTION
template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::MosekRed::MosekRed(const size_t maxLogTerms)
	: MosekBase(NCDim+CDim, CDim, NCDim, NumConstraints+1, maxLogTerms, "red")
{
	mskcall(MSK_putintparam(task, MSK_IPAR_INTPNT_MAX_ITERATIONS, 1e5));
}
#endif


template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::MosekBnd::MosekBnd(const size_t maxLogTerms)
	: MosekBase(NCDim+CDim+1, CDim, NCDim, NumConstraints+1, maxLogTerms, "bnd"),
		vOffsetT(NCDim + CDim)
{
	mskcall(MSK_putvarbound(task, vOffsetT, MSK_BK_FR, -MSK_INFINITY, +MSK_INFINITY));
	MosekBase::putobj(vOffsetT, 1.0);
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
size_t
RR<NCDim,CDim,NumConstraints>::MosekBase::setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign)
{
	if (++logIdx >= maxLog)
		throw std::logic_error(std::string(__func__) + ": too many Log terms");

	// add aux constraint and variable
	mskcall(MSK_appendvars(task, 1));
	mskcall(MSK_appendcons(task, 1));
	putvarbound(numVars+logIdx, 0.0);
	//mskcall(MSK_putvarbound(task, dim+logIdx, MSK_BK_LO, 0.0, +MSK_INFINITY)); // TODO _FR

	// set aux constraint
	assert(linC.size() == dimC);
	assert(dimNC == 0 || linNC.size() == dimNC);

	for (size_t i = 0; i < dimC; ++i)
	{
		if (linC[i] != 0.0)
			putaij(numConstr+logIdx, vOffsetC + i, linC[i]);
	}

	for (size_t i = 0; i < dimNC; ++i)
	{
		if (linNC[i] != 0.0)
			putaij(numConstr+logIdx, vOffsetNC + i, linNC[i]);
	}
	
	putaij(numConstr+logIdx, numVars+logIdx, -1.0);

#if DEBUG > 4
	std::cout << name << ": putconbound(" << numConstr+logIdx << ", 0.0) / Fixed" << std::endl;
#endif
	mskcall(MSK_putconbound(task, numConstr+logIdx, MSK_BK_FX, 0.0, 0.0));

#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << ": setLogConstr(opric = " << constrIdx << ", oprjc = " << numVars+logIdx << ") : " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

	// non-linear term
	update_MSKsch = true;
	oprc[logIdx] = MSK_OPR_LOG;
	opric[logIdx] = constrIdx;
	oprjc[logIdx] = numVars+logIdx;
	oprfc[logIdx] = negativeSign ? 1.0 : -1.0;
	oprgc[logIdx] = 1.0;
	oprhc[logIdx] = c;

	return logIdx;
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekBase::updateModel()
{
	if (update_MSKsch)
	{
#if DEBUG > 4
	const int BT = 3;
	void *buffer[BT];
	char **strings;
	int nptr = backtrace(buffer, BT);
	strings = backtrace_symbols(buffer, nptr);
	std::cout << name << "updateModel : " << strings[2] << " : " << strings[1] << " : " << strings[0] << std::endl;
	free(strings);
#endif

		if (MSKsch != nullptr)
			MSK_scend(task, &MSKsch);

		mskcall(MSK_scbegin(task, 0, nullptr, nullptr, nullptr, nullptr, nullptr, logIdx+1, oprc.get(), opric.get(), oprjc.get(), oprfc.get(), oprgc.get(), oprhc.get(), &MSKsch));

		update_MSKsch = false;
	}
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
RR<NCDim,CDim,NumConstraints>::MosekObj::MosekObj(const size_t maxLogTerms)
	: MosekBase(CDim, CDim, 0, NumConstraints, maxLogTerms, "obj")
{
	mskcall(MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE));
}
	
template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekObj::updateLogConstant(const size_t logIdx, const basetype c)
{
	update_MSKsch = true;
	oprhc[logIdx] = c;
}

template <size_t NCDim, size_t CDim, size_t NumConstraints>
void
RR<NCDim,CDim,NumConstraints>::MosekObj::putobj(const vtypeC& numerLinC, const vtypeC& denomLinC)
{
	numer = numerLinC;
	denom = denomLinC;

	if (!std::all_of(denomLinC.begin(), denomLinC.end(), [] (basetype e) { return e == 0.0; }))
	{
		dinkelbach = true;
	}
	else
	{
		dinkelbach = false;

		for (size_t i = 0; i < CDim; ++i)
			MosekBase::putobj(vOffsetC + i, numerLinC[i]);
	}
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
inline void
RR<NCDim,CDim,NumConstraints>::MosekObj::putobjConst(const basetype numer, const basetype denom)
{
	numerConst = numer;
	denomConst = denom;
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
basetype
RR<NCDim,CDim,NumConstraints>::MosekObj::getobj() const
{
	if (dinkelbach)
	{
		return dinkelbach_lambda;
	}
	else
	{
		basetype sol = MosekBase::getobj();

		return (sol + numerConst) / denomConst;
	}
}


template <size_t NCDim, size_t CDim, size_t NumConstraints>
bool
RR<NCDim,CDim,NumConstraints>::MosekObj::solve()
{
	if (dinkelbach)
	{
		dinkelbach_lambda = 0.0;
		basetype aux = 0.0;

		do
		{
			vtypeC x;

			for (size_t i = 0; i < CDim; ++i)
				MosekBase::putobj(vOffsetC + i, numer[i] - dinkelbach_lambda*denom[i]);

			if ( ! MosekBase::solve())
				return false;

			aux = MosekBase::getobj() + numerConst - dinkelbach_lambda*denomConst;
			getopt_C(x);

			dinkelbach_lambda = std::inner_product(numer.begin(), numer.end(), x.begin(), numerConst) / std::inner_product(denom.begin(), denom.end(), x.begin(), denomConst);
		} while (std::abs(aux) < dinkelbach_tol);

		return true;
	}
	else
		return MosekBase::solve();
}
