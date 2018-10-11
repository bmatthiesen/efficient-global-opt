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

#include "HK.h"

const unsigned short HK::mx_l[] = {2,0,1};
const unsigned short HK::mx_q[] = {1,2,0};


HK::HK(const BndSolvers bs)
	: _RR(dim*10+3, bs, ObjSolvers::Mosek)
{
	this->outputEvery = 10000;
}


void
HK::init(H5::WPPtr wp)
{
	if (initialized)
		throw std::runtime_error("already initialized. create new object instead.");
	initialized = true;

	auto P0 = std::pow(10, wp->P0/10.0);
	decltype(wp->P) P;
	std::transform(wp->P.begin(), wp->P.end(), P.begin(), [](const double p) { return std::pow(10, p/10.0);});

	std::array<basetype, HKDIM> habs;
	std::array<basetype, HKDIM> gabs;

	std::transform(wp->h.begin(), wp->h.end(), habs.begin(), [&wp](const H5::Complex& a) { return abssquare(a) / wp->N0; });
	std::transform(wp->g.begin(), wp->g.end(), gabs.begin(), [&wp,P0](const H5::Complex& a) { return abssquare(a) * P0 / wp->N; });

	// var bounds
	for (size_t i = 0; i < dim; ++i)
	{
		this->setVarbounds_NC(i, 0, P[i]);
		this->setVarbounds_C(i+dim, 0, P[i]);
		this->setVarbounds_C(i, 0);
	}
	this->setVarbounds_NC(dim, 0, std::inner_product(P.begin(), P.end(), habs.begin(), 0.0));

	// some precomputation for logterms
	vtype denom[dim] = {};
	basetype denomConsts[dim] = {};
	vtype A[dim] = {};
	vtype B[dim] = {};
	vtype C[dim] = {};
	vtype D[dim] = {};
	vtypeC AC[dim] = {};
	vtypeC BC[dim] = {};
	vtypeC CC[dim] = {};
	vtypeC DC[dim] = {};

	for (size_t i = 0; i < dim; ++i)
	{
		basetype ginv = 1.0 / gabs[mx_q[i]];

		for (size_t k = 0; k < dim; ++k)
			denom[i][k] = ginv * habs[k];

		denom[i][dim] = ginv;

		denom[i][mx_l[i]] += habs[mx_l[i]];

		denomConsts[i] = 1 + ginv;

		A[i] = denom[i];
		A[i][i] += habs[i];

		B[i] = denom[i];
		B[i][i] += habs[i];
		BC[i][i+dim] = habs[i];

		C[i] = denom[i];
		C[i][i] += habs[i];
		CC[i][mx_l[i]+dim] = habs[mx_l[i]];

		D[i] = denom[i];
		D[i][i] += habs[i];
		DC[i][i+dim] = habs[i];
		DC[i][mx_l[i]+dim] = habs[mx_l[i]];
	}

	// set constraints
	for (size_t i = 0; i < dim; ++i)
	{
		vtypeC ar = {};

		// R_k < B_k
		ar[i] = 1.0;
		this->setConstraint(i, ar, {});
		this->setLogterm(i, BC[i], B[i], denom[i], denomConsts[i], denomConsts[i]);

		// R_k + R_q < A_k + D_q
		ar = {};
		ar[i] = 1.0;
		ar[mx_q[i]] = 1.0;
		setConstraint(dim+i, ar, {});

		this->setLogterm(dim+i, AC[i], A[i], denom[i], denomConsts[i], denomConsts[i]);
		this->setLogterm(dim+i, DC[mx_q[i]], D[mx_q[i]], denom[mx_q[i]], denomConsts[mx_q[i]], denomConsts[mx_q[i]]);

		// R_k + R_q + R_l < A_k + C_q + D_l
		ar = {};
		ar[i] = 1.0;
		ar[mx_q[i]] = 1.0;
		ar[mx_l[i]] = 1.0;
		this->setConstraint(2*dim+i, ar, {});

		this->setLogterm(2*dim+i, AC[i], A[i], denom[i], denomConsts[i], denomConsts[i]);
		this->setLogterm(2*dim+i, CC[mx_q[i]], C[mx_q[i]], denom[mx_q[i]], denomConsts[mx_q[i]], denomConsts[mx_q[i]]);
		this->setLogterm(2*dim+i, DC[mx_l[i]], D[mx_l[i]], denom[mx_l[i]], denomConsts[mx_l[i]], denomConsts[mx_l[i]]);

		// 2*R_k + R_q + R_l < A_k + C_q + C_l + D_k
		ar = {};
		ar[i] = 2.0;
		ar[mx_q[i]] = 1.0;
		ar[mx_l[i]] = 1.0;
		this->setConstraint(3*dim+i, ar, {});

		this->setLogterm(3*dim+i, AC[i], A[i], denom[i], denomConsts[i], denomConsts[i]);
		this->setLogterm(3*dim+i, CC[mx_q[i]], C[mx_q[i]], denom[mx_q[i]], denomConsts[mx_q[i]], denomConsts[mx_q[i]]);
		this->setLogterm(3*dim+i, CC[mx_l[i]], C[mx_l[i]], denom[mx_l[i]], denomConsts[mx_l[i]], denomConsts[mx_l[i]]);
		this->setLogterm(3*dim+i, DC[i], D[i], denom[i], denomConsts[i], denomConsts[i]);

		// pk^c + pk^s <= Pk TODO
		ar = {};
		ar[dim+i] = 1.0;
		vtype br = {};
		br[i] = 1.0;
		this->setConstraint(4*dim+i, ar, br, P[i]);
	}

	{
		// R_1 + R_2 + R_3 < C_1 + C_2 + C_3
		vtypeC ar = {};

		ar[0] = 1.0;
		ar[1] = 1.0;
		ar[2] = 1.0;

		this->setConstraint(5*dim, ar, {});

		this->setLogterm(5*dim, CC[0], C[0], denom[0], denomConsts[0], denomConsts[0]);
		this->setLogterm(5*dim, CC[1], C[1], denom[1], denomConsts[1], denomConsts[1]);
		this->setLogterm(5*dim, CC[2], C[2], denom[2], denomConsts[2], denomConsts[2]);

		// aux
		ar = {};
		std::copy(habs.begin(), habs.end(), ar.begin()+HKDIM);
		vtype br = {};
		br[dim] = -1.0;

		this->setConstraint(5*dim+1, ar, br, 0.0, Boundtype::Upper);
	}
}


void
HK::optimize()
{
	if (!initialized)
		throw std::runtime_error("initialize first");

	_RR::optimize();

	std::cout << "Peak RSS: " << getPeakRSS() << std::endl;
}
