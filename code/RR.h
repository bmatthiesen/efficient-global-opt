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

#ifndef _RR_H
#define _RR_H

#include <limits>
#include <array>
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>

//#include "MosekSep.h"
#include "gurobi_c++.h"

#include "SIT.h"

#define ERR(s)  std::string(__func__) + ": " + (s)

#include "Mosek.h"
#include "Gurobi.h"

enum class BndSolvers { Mosek, MosekSep, GurobiSep };
enum class ObjSolvers { Mosek, Gurobi };

template <size_t NCDim, size_t CDim, size_t NumConstraints>
class RR : protected SIT<NCDim>
{
	public:
		using typename SIT<NCDim>::vtype;
		using vtypeC = std::array<basetype, CDim>;
		using typename SIT<NCDim>::RBox;
		using typename SIT<NCDim>::PBox;

		enum class Boundtype { Undefined = 0, Upper, Fixed, Lower };

		struct Logterm
		{
			const vtypeC linC;
			const vtype linNC;
			const vtype dcLinNC;
			const basetype c;
			const basetype dcC;
			const bool negativeSign;

			size_t obj_logIdx;

			bool allzero_linC;
			bool allzero_linNC;
			bool allzero_dcLinNC;

			Logterm(const vtypeC& _linC, const vtype& _linNC, const vtype& _dcLinNC, const basetype _c, const basetype _dcC, const bool negativeSign)
				: linC(_linC), linNC(_linNC), dcLinNC(_dcLinNC), c(_c), dcC(_dcC), negativeSign(negativeSign),
				allzero_linC(std::all_of(_linC.begin(), _linC.end(), [](basetype e) { return e == 0.0; })),
				allzero_linNC(std::all_of(_linNC.begin(), _linNC.end(), [](basetype e) { return e == 0.0; })),
				allzero_dcLinNC(std::all_of(_dcLinNC.begin(), _dcLinNC.end(), [](basetype e) { return e == 0.0; }))
			{ }
		};

#include "bits/Solver.h"
#include "bits/Gurobi.h"
#include "bits/Mosek.h"

		/* ctor */
		RR(
				const size_t maxLogterms,
				const BndSolvers bndS = BndSolvers::Mosek,
				const ObjSolvers objS = ObjSolvers::Gurobi
		   );

		virtual ~RR() {};

		/* types etc */

		static_assert(std::numeric_limits<basetype>::has_quiet_NaN, "basetype has no NaN");
		static constexpr basetype Unconstrained = std::numeric_limits<basetype>::quiet_NaN();

		/* variables */
		const vtypeC& optC = _optC;
		const vtype& optNC = xopt;

		/* interface */
		void setVarbounds_NC(const size_t varIdx, const basetype lb, const basetype ub);
		void setVarbounds_C (const size_t varIdx, const basetype lb, const basetype ub = Unconstrained);

		void setObjective(const vtypeC& numerLinC, const vtype& numerLinNC, const basetype numerConst = 0.0);
		void setObjective(const vtypeC& numerLinC, const vtype& numerLinNC, const vtypeC& denomLinC, const vtype& denomLinNC, const basetype numerConst = 0.0, const basetype denomConst = 0.0);

		void setConstraint(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c = 0.0, const Boundtype bnd = Boundtype::Upper);
		void setLogterm(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const vtype& dcLinNC, const basetype c = 0.0, const basetype dcC = 0.0, const bool negativeSign = false);

		/* SIT */
		void optimize() override;
		basetype obj(const vtype& P, bool update) override final
			{ return obj(P, nullptr, update); }
		void printResult() const override;

		using SIT<NCDim>::setPrecision;
		using SIT<NCDim>::output;
#ifdef REDUCTION
		using SIT<NCDim>::disableReduction;
#endif

		using SIT<NCDim>::enableBackup;
		using SIT<NCDim>::disableBackup;
		using SIT<NCDim>::removeBackup;

		using SIT<NCDim>::xopt;
		using SIT<NCDim>::optval;
		using SIT<NCDim>::iter;
		using SIT<NCDim>::lastUpdate;
		using SIT<NCDim>::status;
		using SIT<NCDim>::statusStr;
		using SIT<NCDim>::runtime;

		using SIT<NCDim>::Status;

	protected:
		using SIT<NCDim>::gamma;
		using SIT<NCDim>::gamma0;
		using SIT<NCDim>::epsilon;
		using SIT<NCDim>::eta;
		using SIT<NCDim>::setStatus;

		/* variables */
		vtypeC _optC;

		constexpr size_t dimNC() const { return NCDim; }
		constexpr size_t dimC() const { return CDim; }
		constexpr size_t numConstr() const { return NumConstraints; }

		const size_t maxLog;
		const size_t cOffsetObj = 0;
		const size_t cOffsetC = 1;

		std::array<basetype, NumConstraints> constrConsts = {};
		std::array<std::vector<Logterm>, NumConstraints> logterms = {};
		vtype obj_numerLinNC = {};
		vtypeC obj_numerLinC = {};
		vtype obj_denomLinNC = {};
		vtypeC obj_denomLinC = {};
		basetype obj_numerConst = 0.0;
		basetype obj_denomConst = 0.0;
		bool obj_denomPresent = false;
		std::array<vtype, NumConstraints> constr_linNC = {};

#ifdef REDUCTION
		std::unique_ptr<RedSolver> red;
#endif
		std::unique_ptr<BndSolver> bnd;
		std::unique_ptr<ObjSolver> objs;

		bool disablePrintResult = false;

		basetype obj(const vtype& P, vtypeC* R)
			{ return obj(P, R, true); }
		basetype obj(const vtype& P, vtypeC* R, bool update);
		void setObjective();

		/* SIT */
#ifdef REDUCTION
		void reduction(const PBox& P, RBox& red) override;
#endif
		void bound(RBox& red) override final;
		bool isFeasible(const vtype& P) override final; // if you remove final: make sure that obj calls correct isFeasible!
		void setGamma(const basetype g) override final;

	private:
		void checkModel();
		void updateObj(const vtype& P);
};


template <size_t NCDim, size_t CDim, size_t NumConstraints>
std::ostream& operator<< (std::ostream &out, const typename RR<NCDim,CDim,NumConstraints>::Boundtype &point);

#include "bits/RR.cpp"
#include "bits/GurobiBnd.cpp"
#include "bits/GurobiObj.cpp"
#include "bits/Mosek.cpp"
#include "bits/MosekSep.cpp"

#endif
