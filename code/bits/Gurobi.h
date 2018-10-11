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

struct LogtermGrb
{
	const vtype linNC;
	const basetype c;
	const bool negativeSign;

	size_t constrIdx;

	LogtermGrb(const vtype& _linNC, const basetype _c, const size_t constrIdx, const bool negativeSign)
		: linNC(_linNC), c(_c), negativeSign(negativeSign), constrIdx(constrIdx)
	{ }
};

class GurobiBnd : public BndSolver
{
	private:
		GRBModel Grb;
		std::array<GRBConstr, NumConstraints+1> constr_linC = {};
		std::unique_ptr<GRBVar[]> varsC;
		GRBVar varT;

		std::array<basetype, NumConstraints+1> constr_consts = {};
		std::array<vtype, NumConstraints+1> constr_linNC = {};
		std::vector<LogtermGrb> logterms;
		vtype NCub = {};
		std::array<short int, NCDim> NCubSgnU = {};

	public:
		GurobiBnd(const size_t maxLogTerms);

		void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub)
			{ varsC[varIdx].set(GRB_DoubleAttr_LB, lb); varsC[varIdx].set(GRB_DoubleAttr_UB, ub); }

		void putvarbound_C(const size_t varIdx, const basetype lb)
			{ varsC[varIdx].set(GRB_DoubleAttr_LB, lb); }

		void putvarbound_NC(const size_t varIdx, const basetype lb, const basetype ub);

		void putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val);

		void putaij_C(const size_t constrIdx, const size_t varIdx, const basetype val)
			{ Grb.chgCoeff(constr_linC[constrIdx], varsC[varIdx], val); }

		void putaij_T(const size_t constrIdx, const basetype val)
			{ Grb.chgCoeff(constr_linC[constrIdx], varT, val); }

		void putconbound(const size_t constrIdx, const basetype ub)
			{ constr_consts[constrIdx] = ub; }

		basetype getobj() const
			{ return Grb.get(GRB_DoubleAttr_ObjVal); }

		bool solve();

		void setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd);
		size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign);

		void updateModel();

		/* BndSolver */
		void getopt_NC(vtype& x) const
			{ x = NCub; }
};

class GurobiObj : public ObjSolver
{
	private:
		GRBModel Grb;
		std::array<GRBConstr, NumConstraints> constr_linC = {};
		std::unique_ptr<GRBVar[]> varsC;

		basetype numerConst = 0;
		basetype denomConst = 0;


	public:
		GurobiObj();

		void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub)
			{ varsC[varIdx].set(GRB_DoubleAttr_LB, lb); varsC[varIdx].set(GRB_DoubleAttr_UB, ub); }

		void putvarbound_C(const size_t varIdx, const basetype lb)
			{ varsC[varIdx].set(GRB_DoubleAttr_LB, lb); }

		void putconbound(const size_t constrIdx, const basetype ub)
			{ constr_linC[constrIdx].set(GRB_DoubleAttr_RHS, ub); }

		void setConstr(const size_t constrIdx, const vtypeC& linC, const basetype c, const Boundtype bnd);
		size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const basetype c, const bool negativeSign);
		void updateLogConstant(const size_t logIdx, const basetype c);

		void putobj(const vtypeC& numerLinC, const vtypeC& denomLinC);
		void putobjConst(const basetype numer, const basetype denom);

		basetype getobj() const;

		void updateModel()
			{ Grb.update(); }

		bool solve();

		void getopt_C(vtypeC& x) const;
};
