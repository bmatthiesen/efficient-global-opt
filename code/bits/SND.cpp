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

template <unsigned long long sndtype>
const unsigned short SND<sndtype>::mx_l[] = {2,0,1};
template <unsigned long long sndtype>
const unsigned short SND<sndtype>::mx_q[] = {1,2,0};


template <unsigned long long sndtype>
SND<sndtype>::SND(const BndSolvers bs)
	: _RR(dim*2, bs, ObjSolvers::Gurobi), useSND(sndtype)
{
	_RR::outputEvery = 10000;
}


template <unsigned long long sndtype>
void
SND<sndtype>::init(H5::WPPtr wp, const bool useTighterBound)
{
	if (initialized)
		throw std::runtime_error("already initialized. create new object instead.");
	initialized = true;


	auto P0 = std::pow(10, wp->P0/10.0);
	decltype(wp->P) P;
	std::transform(wp->P.begin(), wp->P.end(), P.begin(), [](const double p) { return std::pow(10, p/10.0);});

	H5::vec_t habs;
	H5::vec_t gabs;

	std::transform(wp->h.begin(), wp->h.end(), habs.begin(), [&wp](const H5::Complex& a) { return abssquare(a) / wp->N0; });
	std::transform(wp->g.begin(), wp->g.end(), gabs.begin(), [&wp,P0](const H5::Complex& a) { return abssquare(a) * P0 / wp->N; });

	for (size_t i = 0; i < dim; ++i)
	{
		_RR::setVarbounds_NC(i, 0, P[i]);
		_RR::setVarbounds_C(i, 0);
	}

	size_t cidx = -1;
	for (size_t i = 0; i < dim; ++i)
	{
		vtypeC ar {}; // {} initializes to zero
		vtype f {};
		vtype g {};

		basetype ginv = 1.0 / gabs[mx_q[i]];

		std::transform(habs.begin(), habs.end(), g.begin(), [&ginv] (const basetype& h) { return ginv * h; });
		ginv += 1;

		if (useSND[i])
		{
			 /* single user */
			ar[i] = 1.0;
			f = g;
			f[i] += habs[i];

			this->setConstraint(++cidx, ar, {});

			if (useTighterBound)
				this->setLogterm(cidx, {}, f, g, ginv, ginv, false);
			else
				this->setLogterm(cidx, {}, g, f, ginv, ginv, true);

			/* multi user */
			ar[mx_l[i]] = 1.0;
			f[mx_l[i]] += habs[mx_l[i]];

			this->setConstraint(++cidx, ar, {});

			if (useTighterBound)
				this->setLogterm(cidx, {}, f, g, ginv, ginv, false);
			else
				this->setLogterm(cidx, {}, g, f, ginv, ginv, true);
		}
		else
		{
			g[mx_l[i]] += habs[mx_l[i]];
			ar[i] = 1.0;
			f = g;
			f[i] += habs[i];

			this->setConstraint(++cidx, ar, {});

			if (useTighterBound)
				this->setLogterm(cidx, {}, f, g, ginv, ginv, false);
			else
				this->setLogterm(cidx, {}, g, f, ginv, ginv, true);
		}
	}
}

template <unsigned long long sndtype>
void
SND<sndtype>::optimize()
{
	if (!initialized)
		throw std::runtime_error("initialize first");

	_RR::optimize();

	std::cout << "Peak RSS: " << getPeakRSS() << std::endl;
}

