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

template <size_t numUE, size_t numNC, bool simple>
double SNDBench<numUE,numNC,simple>::N = 1e-2;

template <size_t numUE, size_t numNC, bool simple>
double SNDBench<numUE,numNC,simple>::P = 1;

template <size_t numUE, size_t numNC, bool simple>
SNDBench<numUE, numNC,simple>::SNDBench(const BndSolvers bs)
	: _RR(numConstr, bs, ObjSolvers::Gurobi)
{
	this->outputEvery = 10000;
}

template <size_t numUE, size_t numNC, bool simple>
void
SNDBench<numUE, numNC,simple>::init(typename H5::WPPtr wp)
{
	if (initialized)
		throw std::runtime_error("already initialized. create new object instead.");
	initialized = true;

	/* effective channel */
	std::array<basetype, numUE> h[numUE];

	for (size_t i = 0; i < numUE; ++i)
		for (size_t j = 0; j < numUE; ++j)
			h[i][j] = abssquare(wp->h[i][j]);

	/* set variable bounds */
	for (size_t i = 0; i < numUE; ++i)
	{
		this->setVarbounds_C(oR + i, 0); // rates

		if (i < numNC)
			this->setVarbounds_NC(i, 0, P);
	}
	
	/* create base sets for rate region */
	std::set<size_t> cSet, ncSet;
	for (size_t i = 0; i < numNC; ++i)
		ncSet.insert(i);
	for (size_t i = numNC; i < numUE; ++i)
		cSet.insert(i);

	/* rate regions */
	size_t cidx = -1; // constraint idx
	for (size_t i = 0; i < numUE; ++i)
	{
#if DEBUG > 2
		std::cout << "RR for " << i << std::endl;
#endif

		/* create sets S and Sc */
		std::set<size_t> s;

		if constexpr(simple)
			s = {i, (i+1)%numUE};
		else
			s = cSet;

		std::set<size_t> sc(ncSet);

		if constexpr(simple)
		{
			for (auto&& j: s)
			{
				decltype(sc)::const_iterator it;
				if ((it = sc.find(j)) != sc.end())
					sc.erase(it);
			}
		}
		else
		{
			if (i < numNC)
			{
				sc.erase(i);
				s.insert(i);
			}
		}

#if DEBUG > 2
		std::cout << "\tS = " << s << " ; Sc = " << sc << std::endl;
#endif

		/* P_k(Sc) = denom */
		vtype denom = {};
		for (auto&& j: sc)
			denom[j] = h[i][j];

		/* create power sets T */
		auto pset = powerset(s);

		for (auto&& t: pset)
		{
			if (t.size() == 0)
				continue;

			vtypeC a = {};
			double numerC = 0;
			vtype numerNC = denom;

			for (auto&& u: t)
			{
				a[oR + u] = 1;

				if (u < numNC)
					numerNC[u] += h[i][u];
				else
					// convex variable: index is oP - numNC + u
					numerC += h[i][u] * P;
			}

			this->setConstraint(++cidx, a, {});
			this->setLogterm(cidx, {}, numerNC, denom, N+numerC, N, false);

#if DEBUG > 2
			std::cout << "\t\t" << cidx << ": " << t << std::endl;
#endif
		}

#if DEBUG > 2
		std::cout << std::endl;
#endif
	}
}

template <size_t numUE, size_t numNC, bool simple>
void
SNDBench<numUE, numNC,simple>::optimize()
{
	if (!initialized)
		throw std::runtime_error("initialize first");

	_RR::optimize();

	std::cout << "Peak RSS: " << getPeakRSS() << std::endl;
}

// source: https://rosettacode.org/mw/index.php?title=Power_set&oldid=278114#C.2B.2B14_version
template <class S>
auto powerset(const S& s)
{
    std::set<S> ret;
    ret.emplace();
    for (auto&& e: s)
	{
        std::set<S> rs;
        for (auto x: ret)
		{
            x.insert(e);
            rs.insert(x);
        }
        ret.insert(begin(rs), end(rs));
    }
    return ret;
}

// set printer
template <typename T>
std::ostream& operator<<(std::ostream &out, const std::set<T>& s)
{
	out << "{ ";
	char const* prefix = "";
	for (auto&& e: s) {
		out << prefix << e;
		prefix = ", ";
	}
	out << " }";

	return out;
}
