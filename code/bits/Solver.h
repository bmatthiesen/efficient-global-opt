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

class BndSolver
{
	public:
		virtual ~BndSolver() noexcept {};

		virtual void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub) =0;
		virtual void putvarbound_C(const size_t varIdx, const basetype lb) =0;
		virtual void putvarbound_NC(const size_t varIdx, const basetype lb, const basetype ub) =0;

		virtual void putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val) =0;
		virtual void putaij_C(const size_t constrIdx, const size_t varIdx, const basetype val) =0;


		virtual void putconbound(const size_t constrIdx, const basetype ub) =0;
		virtual basetype getobj() const =0;
		virtual bool solve() =0;

		virtual void setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd) =0;

		virtual size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign) =0;
		virtual void updateModel() =0;

		virtual void getopt_NC(vtype& x) const =0;
		virtual void putaij_T(const size_t constrIdx, const basetype val) =0;
};

#ifdef REDUCTION
class RedSolver
{
	public:
		virtual ~RedSolver() noexcept {};

		virtual void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub) =0;
		virtual void putvarbound_C(const size_t varIdx, const basetype lb) =0;
		virtual void putvarbound_NC(const size_t varIdx, const basetype lb, const basetype ub) =0;

		virtual void putaij_NC(const size_t constrIdx, const size_t varIdx, const basetype val) =0;
		virtual void putaij_C(const size_t constrIdx, const size_t varIdx, const basetype val) =0;


		virtual void putconbound(const size_t constrIdx, const basetype ub) =0;
		virtual basetype getobj() const =0;
		virtual bool solve() =0;

		virtual void setConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const Boundtype bnd) =0;

		virtual size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const vtype& linNC, const basetype c, const bool negativeSign) =0;
		virtual void updateModel() =0;

		virtual void setMinimize() =0;
		virtual void setMaximize() =0;

		virtual void putobj_NC(const size_t varIdx, const basetype val) =0;
};
#endif


class ObjSolver
{
	public:
		virtual ~ObjSolver() noexcept {};

		virtual void putvarbound_C(const size_t varIdx, const basetype lb, const basetype ub) =0;
		virtual void putvarbound_C(const size_t varIdx, const basetype lb) =0;

		virtual void putconbound(const size_t constrIdx, const basetype ub) =0;
		virtual void setConstr(const size_t constrIdx, const vtypeC& linC, const basetype c, const Boundtype bnd) =0;
		virtual size_t setLogConstr(const size_t constrIdx, const vtypeC& linC, const basetype c, const bool negativeSign) =0;
		virtual void updateLogConstant(const size_t logIdx, const basetype c) =0;

		virtual void putobj(const vtypeC& numerLinC, const vtypeC& denomLinC) =0;
		virtual void putobjConst(const basetype numer, const basetype denom) =0;
		virtual basetype getobj() const =0;

		virtual void updateModel() =0;
		virtual bool solve() =0;

		virtual void getopt_C(vtypeC& x) const =0;
};
