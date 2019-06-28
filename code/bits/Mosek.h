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

// included from RR.h

struct LogtermMskbnd
{
	const vtype linNC;
	const basetype c;

	size_t constrIdx;

	LogtermMskbnd(const vtype& _linNC, const basetype _c, const size_t constrIdx)
		: linNC(_linNC), c(_c), constrIdx(constrIdx)
	{ }
};


class MosekBase
{
	protected:
		MSKtask_t task;
		std::vector<Boundtype> btypes;

		const size_t numVars;
		const size_t dimNC;
		const size_t dimC;
		size_t logIdx;
		const size_t numConstr;
		const size_t maxLog;

		const size_t vOffsetNC;
		const size_t vOffsetC;

		schand_t MSKsch;
		bool update_MSKsch;

		std::unique_ptr<int[]> oprc;
		std::unique_ptr<int[]> opric;
		std::unique_ptr<int[]> oprjc;

		std::unique_ptr<double[]> oprfc;
		std::unique_ptr<double[]> oprgc;
		std::unique_ptr<double[]> oprhc;

		const char* name;


		MosekBase(const size_t numVars, const size_t Cdim, const size_t NCdim, const size_t numConstraints, const size_t maxLogTerms, const char* name);
		virtual ~MosekBase() noexcept;

		void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub)
			{ putvarbound(vOffsetC + varIdx, lb, ub); }

		void putvarbound_C(const size_t varIdx, const basetype lb)
			{ putvarbound(vOffsetC + varIdx, lb); }

		void putvarbound_NC(const size_t varIdx, const basetype lb, const basetype ub)
			{ putvarbound(vOffsetNC + varIdx, lb, ub); }

		void putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val)
			{ putaij(constrIdx, vOffsetNC + varIdx, val); }

		void putaij_C(const size_t constrIdx, const size_t varIdx, const basetype val)
			{ putaij(constrIdx, vOffsetC + varIdx, val); }

		void putconbound(const size_t constrIdx, const basetype ub);
		basetype getobj() const;
		bool solve();

		void setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd);

		size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign);
		void updateModel();

		void putvarbound(const size_t varIdx, const basetype lb, const basetype ub);
		void putvarbound(const size_t varIdx, const basetype lb);
		void putaij(const size_t constrIdx, const size_t varIdx, const basetype val);
		void putobj(const size_t varIdx, const basetype val);
		void getopt(basetype* x, const size_t offset, const size_t len) const;

	private:
		void destroy() noexcept;
};

#ifdef REDUCTION
class MosekRed : private MosekBase, public RedSolver
{
	protected:
		using MosekBase::task;
		using MosekBase::vOffsetNC;
		using MosekBase::vOffsetC;

	public:
		MosekRed(const size_t maxLogTerms);

		/* SolverBase */
		void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub)
			{ MosekBase::putvarbound_C(varIdx, lb, ub); }

		void putvarbound_C(const size_t varIdx, const basetype lb)
			{ MosekBase::putvarbound_C(varIdx, lb); }

		void putvarbound_NC(const size_t varIdx, const basetype lb, const basetype ub)
			{ MosekBase::putvarbound_NC(varIdx, lb, ub); }

		void putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val)
			{ MosekBase::putaij_NC(constrIdx, varIdx, val); }

		void putaij_C(const size_t constrIdx, const size_t varIdx, const basetype val)
			{ MosekBase::putaij_C(constrIdx, varIdx, val); }

		void putconbound(const size_t constrIdx, const basetype ub)
			{ MosekBase::putconbound(constrIdx, ub); }

		basetype getobj() const
			{ return MosekBase::getobj(); }

		bool solve()
			{ return MosekBase::solve(); }

		void setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd)
			{ MosekBase::setConstr(constrIdx, linC, linNC, c, bnd); }

		size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign)
			{ return MosekBase::setLogConstr(constrIdx, linC, linNC, c, negativeSign); }

		void updateModel()
			{ MosekBase::updateModel(); }

		/* RedSolver */
		void setMinimize();
		void setMaximize();

		void putobj_NC(const size_t varIdx, const basetype val)
			{ MosekBase::putobj(vOffsetNC + varIdx, val); }
};
#endif

class MosekBnd : private MosekBase, public BndSolver
{
	protected:
		const size_t vOffsetT;

		using MosekBase::task;
		using MosekBase::vOffsetNC;
		using MosekBase::vOffsetC;

	public:
		MosekBnd(const size_t maxLogTerms);

		/* SolverBase */
		void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub)
			{ MosekBase::putvarbound_C(varIdx, lb, ub); }

		void putvarbound_C(const size_t varIdx, const basetype lb)
			{ MosekBase::putvarbound_C(varIdx, lb); }

		void putvarbound_NC(const size_t varIdx, const basetype lb, const basetype ub)
			{ MosekBase::putvarbound_NC(varIdx, lb, ub); }

		void putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val)
			{ MosekBase::putaij_NC(constrIdx, varIdx, val); }

		void putaij_C(const size_t constrIdx, const size_t varIdx, const basetype val)
			{ MosekBase::putaij_C(constrIdx, varIdx, val); }

		void putconbound(const size_t constrIdx, const basetype ub)
			{ MosekBase::putconbound(constrIdx, ub); }

		basetype getobj() const
			{ return MosekBase::getobj(); }

		bool solve()
			{ return MosekBase::solve(); }

		void setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd)
			{ MosekBase::setConstr(constrIdx, linC, linNC, c, bnd); }

		size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign)
			{ return MosekBase::setLogConstr(constrIdx, linC, linNC, c, negativeSign); }

		void updateModel()
			{ MosekBase::updateModel(); }

		/* BndSolver */
		void getopt_NC(vtype& x) const
			{ this->getopt(x.data(), vOffsetNC, NCDim); }

		void putaij_T(const size_t constrIdx, const basetype val)
			{ MosekBase::putaij(constrIdx, vOffsetT, val); }
};

class MosekSepBnd : private MosekBase, public BndSolver
{
	private:
		const size_t vOffsetT;
		std::array<basetype, NumConstraints> constr_consts;
		std::array<vtype, NumConstraints> constr_linNC;
		std::vector<LogtermMskbnd> logterms;
		vtype NCub;

		using MosekBase::task;

	public:
		MosekSepBnd(const size_t maxLogTerms);

		/* SolverBase */
		void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub)
			{ MosekBase::putvarbound_C(varIdx, lb, ub); }

		void putvarbound_C(const size_t varIdx, const basetype lb)
			{ MosekBase::putvarbound_C(varIdx, lb); }

		void putvarbound_NC(const size_t varIdx, const basetype, const basetype ub)
			{ NCub[varIdx] = ub; }

		void putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val);

		void putaij_C(const size_t constrIdx, const size_t varIdx, const basetype val)
			{ MosekBase::putaij_C(constrIdx, varIdx, val); }

		void putconbound(const size_t constrIdx, const basetype ub);

		basetype getobj() const
			{ return MosekBase::getobj(); }

		bool solve()
			{updateModel(); return MosekBase::solve(); }

		void setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd);
		size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign);

		void updateModel();

		/* BndSolver */
		void getopt_NC(vtype& x) const
			{ x = NCub; }

		void putaij_T(const size_t constrIdx, const basetype val)
			{ MosekBase::putaij(constrIdx, vOffsetT, val); }
};

class MosekObj : private MosekBase, public ObjSolver
{
	protected:
		using MosekBase::task;
		using MosekBase::vOffsetNC;
		using MosekBase::vOffsetC;
		using MosekBase::oprhc;
		using MosekBase::update_MSKsch;

	private:
		basetype numerConst = 0;
		basetype denomConst = 0;
		vtypeC numer = {};
		vtypeC denom = {};
		bool dinkelbach = false;
		basetype dinkelbach_lambda;

	public:
		basetype dinkelbach_tol = 1e-8;

		MosekObj(const size_t maxLogTerms);
		void updateLogConstant(const size_t logIdx, const basetype c);

		void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub)
			{ MosekBase::putvarbound_C(varIdx, lb, ub); }

		void putvarbound_C(const size_t varIdx, const basetype lb)
			{ MosekBase::putvarbound_C(varIdx, lb); }

		void putconbound(const size_t constrIdx, const basetype ub)
			{ MosekBase::putconbound(constrIdx, ub); }

		basetype getobj() const;

		bool solve();

		void updateModel()
			{ MosekBase::updateModel(); }

		void setConstr(const size_t constrIdx, const vtypeC& linC, const basetype c, const Boundtype bnd)
			{ MosekBase::setConstr(constrIdx, linC, {}, c, bnd); }

		size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const basetype c, const bool negativeSign)
			{ return MosekBase::setLogConstr(constrIdx, linC, {}, c, negativeSign); }

		void putobj(const vtypeC& numerLinC, const vtypeC& denomLinC);
		void putobjConst(const basetype numer, const basetype denom);

		void getopt_C(vtypeC& x) const
			{ this->getopt(x.data(), vOffsetC, CDim); }
};

