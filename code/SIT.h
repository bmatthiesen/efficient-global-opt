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

#ifndef _SIT_H
#define _SIT_H

#include <memory>
#include <vector>
#include <array>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <stdexcept>
#include <cmath>
#include <queue>
#include <fstream>

#include "SIT.h"
#include "util.h"

using std::cout;
using std::endl;

typedef double basetype;

template <typename T>
std::ostream& operator<< (std::ostream &out, const std::vector<T>& v)
{
	out << "[";
	for (auto& e : v)
		out << " " << e;
	out << " ]";

	return out;
}

template <typename T, size_t N>
std::ostream& operator<< (std::ostream &out, const std::array<T,N>& v)
{
	out << "[";
	for (auto& e : v)
		out << " " << e;
	out << " ]";

	return out;
}


template <typename SIT>
struct compare_RBox
{
	bool operator()(const typename SIT::RBox& n1, const typename SIT::RBox& n2) const
		{ return n1.beta > n2.beta; }
};


template <size_t Dim>
class SIT
{
public:
	using vtype = std::array<basetype, Dim>;

	struct PBox
	{
		vtype lb, ub;
	};

	struct RBox_data
	{
		vtype lb, ub, xk;
		vtype* yk;
	};

	class RBox
	{
		std::unique_ptr<RBox_data> data_;
	public:
		RBox() : data_(std::make_unique<RBox_data>()) {};
		RBox(std::unique_ptr<RBox_data>&& data) : data_(std::move(data)) {};

		std::unique_ptr<RBox_data> move_data() {
			return std::move(data_);
		}

		basetype  lb(size_t index) const { return data_->lb[index]; }
		basetype& lb(size_t index) { return data_->lb[index]; }
		const vtype& lb() const { return data_->lb; }
		vtype& lb() { return data_->lb; }

		basetype  ub(size_t index) const { return data_->ub[index]; }
		basetype& ub(size_t index) { return data_->ub[index]; }
		const vtype& ub() const { return data_->ub; }
		vtype& ub() { return data_->ub; }

		basetype  xk(size_t index) const { return data_->xk[index]; }
		basetype& xk(size_t index) { return data_->xk[index]; }
		const vtype& xk() const { return data_->xk; }
		vtype& xk() { return data_->xk; }

		basetype  yk(size_t index) const { return (*data_->yk)[index]; }
		basetype& yk(size_t index) { return (*data_->yk)[index]; }
		//const vtype& yk() const { return *data_->yk; }
		vtype*& yk() { return data_->yk; }
		const vtype* yk() const { return data_->yk; }

		basetype beta;
	};

	using RType_base = std::priority_queue<RBox, std::vector<RBox>, compare_RBox<SIT>>;

	class RType : public RType_base
	{
		enum class ykptr { undef = 0, xk, lb, ub, null };

	public:
		template <typename... Ts>
		RType(Ts&&... args) : RType_base(std::forward<Ts>(args)...)
		{}
		
		bool backup(std::ofstream& fout);
		bool restore(std::ifstream& fin, MiniPool<RBox_data>& pool);
	};

	class PType
	{
	public:
		PType() : len(1) { };

		PBox& operator[](const size_t index)
		{
			if (len == 1)
				return P1;
			else
				return P2[index];
		}

		const PBox& operator[](const size_t index) const
		{
			if (len == 1)
				return P1;
			else
				return P2[index];
		}

		void use1() { len = 1; }
		void use2() { len = 2; }
		size_t size() { return len; }

	private:
		PBox P1;
		std::array<PBox, 2> P2;
		size_t len;
	};

	enum class Status { Optimal, Unsolved, Infeasible };


	SIT();
	virtual ~SIT() {};

	// parameter setter
	void setPrecision(const basetype eta, const basetype epsilon = -1);
	void setLB(const vtype& v);
	void setLB(const basetype e);
	void setLB(const size_t idx, const basetype e);
	void setUB(const vtype& v);
	void setUB(const basetype e);
	void setUB(const size_t idx, const basetype e);
	void setUserData(const void* ud);

	// parameter getter
	basetype getEta() { return eta; }
	basetype getEpsilon() { return epsilon; }

	// parameter
	basetype gamma0;
	bool output;
	unsigned long long outputEvery;
#ifdef REDUCTION
	bool disableReduction;
#endif

	constexpr size_t dim() const
		{ return Dim; }

	// result
	vtype xopt;
	basetype optval;
	unsigned long long iter, lastUpdate;
	Status status;
	const char *statusStr;
	double runtime; // in seconds

	// run algorithm
	virtual void optimize();

	basetype obj(const vtype& P)
		{ return obj(P, true); }

	virtual basetype obj(const vtype& P, bool update) =0;
	virtual void printResult() const;

	// backup functionality interface
	void enableBackup(const std::string& filename, const std::chrono::seconds& interval);
	void disableBackup()
		{ backupEnable = false; }
	void removeBackup() const;

private:
	// parameter
	vtype lb, ub;

	// backup config
	std::chrono::seconds backupInterval;
	bool backupEnable = false;
	std::string backupFile;
	using backupClock = std::chrono::steady_clock;
	backupClock::time_point backupNext;

	// functions
#ifdef REDUCTION
	virtual void reduction(const PBox& P, RBox& red) =0;
#endif
	virtual void bound(RBox& red) =0;
	
	void initBackup();
	void backup(const std::string& fn, SIT<Dim>::RType& R, SIT<Dim>::PType& P);
	bool restore(const std::string& fn, SIT<Dim>::RType& R, SIT<Dim>::PType& P, MiniPool<RBox_data>& pool);

	basetype maxObj(void);

protected:
	// parameter
	basetype gamma;
	basetype epsilon, eta;
	const void* userdata;

	using clock = std::chrono::high_resolution_clock;
	std::chrono::time_point<clock> tic;

	// functions
	virtual bool isFeasible(const vtype& P) =0;
	virtual void setGamma(const basetype g);
	void setStatus(const Status s);
};


template <size_t Dim>
SIT<Dim>::SIT()
	: gamma0(NAN), output(true), outputEvery(1000),
#ifdef REDUCTION
		disableReduction(true),
#endif
		gamma(0), epsilon(1e-5), eta(1e-2)
{
	setStatus(Status::Unsolved);
}

template <size_t Dim>
void
SIT<Dim>::setPrecision(const basetype eta, const basetype epsilon)
{
	this->eta = eta;
	this->epsilon = epsilon < 0 ? eta * 1e-3 : epsilon;
	setStatus(Status::Unsolved);
}

template <size_t Dim>
void
SIT<Dim>::setLB(const vtype& v)
{
	lb = v;
	setStatus(Status::Unsolved);
}

template <size_t Dim>
void
SIT<Dim>::setLB(const basetype e)
{
	for (auto &b : lb)
		b = e;

	setStatus(Status::Unsolved);
}

template <size_t Dim>
void
SIT<Dim>::setLB(const size_t idx, const basetype e)
{
	lb[idx] = e;
}

template <size_t Dim>
void
SIT<Dim>::setUB(const vtype& v)
{
	ub = v;
	setStatus(Status::Unsolved);
}

template <size_t Dim>
void
SIT<Dim>::setUB(const basetype e)
{
	for (auto &b : ub)
		b = e;

	setStatus(Status::Unsolved);
}

template <size_t Dim>
void
SIT<Dim>::setUB(const size_t idx, const basetype e)
{
	ub[idx] = e;
}

template <size_t Dim>
void
SIT<Dim>::setUserData(const void* ud)
{
	userdata = ud;
	setStatus(Status::Unsolved);
}

template <size_t Dim>
void
SIT<Dim>::printResult() const
{
	cout << "Status: " << statusStr << endl;
	cout << "Optval: " << optval << endl;

	cout << "X*: [";
	std::for_each(xopt.begin(), xopt.end(), [] (const basetype &a) { cout << " " << a; });
	cout << " ]" << endl;

	cout << "Precision: epsilon = " << epsilon << " / eta = " << eta << endl;

	cout << "Iter: " << iter << endl;
	cout << "Solution found in iter: " << lastUpdate << endl;

	cout << "Runtime: " << runtime << " sec" << endl;
}

// SIT Algorithm from Tuy2016, p.209
template <size_t Dim>
void
SIT<Dim>::optimize()
{
	PType P;
	RType R;
	MiniPool<RBox_data> pool;

	runtime = 0.0;

	if (std::isnan(gamma0))
		throw std::runtime_error("Set gamma0!");

	if (backupEnable && restore(backupFile, R, P, pool))
		setGamma(optval + eta);
	else
	{
		// step 0
		iter = lastUpdate = 0;
		optval = -std::numeric_limits<basetype>::infinity();
		setGamma(gamma0);

		P.use1();
		P[0].lb = lb;
		P[0].ub = ub;
	}

	tic = clock::now();
	setStatus(Status::Unsolved);
	initBackup();

	while (true)
	{
		if (output && iter != 0 && iter % outputEvery == 0)
			std::printf("%8llu %8lu %11g  (%llu) | Peak RSS: %zu\n", iter, R.size(), optval, lastUpdate, getPeakRSS());

		if (backupEnable && backupClock::now() >= backupNext)
		{
			backup(backupFile, R, P);
			backupNext = backupClock::now() + backupInterval;
		}

		iter++;

		// step 1: reduce & bound
		for (size_t i = 0; i < P.size(); i++)
		{
			RBox red(pool.get());

#ifdef REDUCTION
			if (disableReduction)
			{
#endif
				red.lb() = std::move(P[i].lb);
				red.ub() = std::move(P[i].ub);
#ifdef REDUCTION
			}
			else
				reduction(P[i], red); // update lb, ub
#endif

			bound(red); // set beta, xk, yk

			if (red.beta <= -epsilon)
			{
#if DEBUG
				bool ret(false);

				for (auto u = red.ub().begin(), l = red.lb().begin(); u != red.ub().end(); ++u, ++l)
					if (*u < *l)
					{
						ret = true;
						break;
					}

				if (ret)
					cout << "this is not supposed to happen" << endl;
#endif

				R.push(std::move(red));
			}
			else
				pool.put(red.move_data());
		}

		// step 2: terminate
		if (R.empty())
		{
			if (gamma == gamma0)
			{
				setStatus(Status::Infeasible);
			}
			else
			{
				optval = obj(xopt); // TODO necessary?
				setStatus(Status::Optimal);
			}

			break;
		}

		// step 3: select box
		{
			typename RType::const_reference M = R.top(); // argmin

			if (isFeasible(M.xk()))
			{
				// step 4: xk is a solution satisfying f(xk) >= gamma
				basetype fxk = obj(M.xk(), false);

				if (gamma == gamma0 || fxk > optval) // TODO test for gamma == gamma0 necessary?
				{
					xopt = M.xk();
					setGamma(fxk + eta);
					optval = std::move(fxk);
					lastUpdate = iter;
				}
			}

			// step 5: branch
			{
				size_t jk = 0;
				basetype max = fabs(M.yk(0) - M.xk(0));
				for (size_t i = 1; i < Dim; i++) // argmax
				{
					basetype diff = fabs(M.yk(i) - M.xk(i));

					if (diff > max)
					{
						jk = i;
						max = std::move(diff);
					}
				}

				basetype vk = (M.yk(jk) + M.xk(jk)) / 2.0;

				P.use2();

				P[0].lb = M.lb();
				P[0].ub = M.ub();
				P[0].ub[jk] = vk;

				P[1].lb = M.lb();
				P[1].lb[jk] = vk;
				P[1].ub = M.ub();

				pool.put(const_cast<RBox&>(M).move_data());
				R.pop();
			}
		}
	}

	runtime += std::chrono::duration<double>(clock::now() - tic).count();

	if (output)
	{
		std::printf("%8llu %8lu %11g  (%llu)\n\n", iter, R.size(), optval, lastUpdate);
		printResult();
	}
}

#if 0
template <size_t Dim>
basetype
SIT<Dim>::maxObj(void)
{
	basetype max = -std::numeric_limits<basetype>::infinity();

	for (int i = 0; i < (1<<Dim); i++)
	{
		vtype c;

		for (size_t j = 0; j < Dim; j++)
			c[j] = i & (1<<j) ? ub[j] : lb[j];

		basetype val = obj(c);

		if (val > max)
			max = val;
	}

	return max;
}
#endif

template <size_t Dim>
void
SIT<Dim>::setGamma(const basetype g)
{
	gamma = g;
}


template <size_t Dim>
void
SIT<Dim>::setStatus(const Status s)
{
	status = s;

	switch (s)
	{
		case Status::Optimal:
			statusStr = "Optimal";
			break;

		case Status::Unsolved:
			statusStr = "Unsolved";
			break;

		case Status::Infeasible:
			statusStr = "Infeasible";
			break;
	}
}


template <size_t Dim>
void
SIT<Dim>::backup(const std::string& fn, SIT<Dim>::RType& R, SIT<Dim>::PType& P)
{
	const auto toc = clock::now();
	bool fail = true;
	const std::string fnbak = fn + ".bak";

	// save runtime since last backup
	runtime += std::chrono::duration<double>(clock::now() - tic).count(); 

	// create backup
	std::remove(fnbak.c_str());
	std::rename(fn.c_str(), fnbak.c_str());

	std::ofstream fout {fn, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc};

	if (!fout)
	{
		std::cerr << "Couldn't open file '" << fn << "' for writing" << std::endl;
		fail = true;
	}
	else
	{
		fout.write(reinterpret_cast<const char*>(&optval), sizeof(optval));
		fout.write(reinterpret_cast<const char*>(&xopt), sizeof(xopt));
		fout.write(reinterpret_cast<const char*>(&runtime), sizeof(runtime));
		fout.write(reinterpret_cast<const char*>(&iter), sizeof(iter));
		fout.write(reinterpret_cast<const char*>(&lastUpdate), sizeof(lastUpdate));

		fout.write(reinterpret_cast<const char*>(&P), sizeof(P));

		if (!R.backup(fout))
			fail = true;
		else
			fail = false;

		fout.close();
	}

	if (fail) // restore backup
	{
		std::remove(fn.c_str());
		std::rename(fnbak.c_str(), fn.c_str());
	}
	else // delete backup
	{
		std::remove(fnbak.c_str());
	}

	if (output && !fail)
	{
		std::cout << ">>> Backup successful: saved " << R.size() << " elements in ";
		std::cout << std::chrono::duration<double>(clock::now() - toc).count() << " seconds" << std::endl;
	}
	else
		std::cerr << ">>> Backup failed" << std::endl;

	tic = clock::now(); // restart tic (this is to exclude backup time from runtime)
}


template <size_t Dim>
bool
SIT<Dim>::restore(const std::string& fn, SIT<Dim>::RType& R, SIT<Dim>::PType& P, MiniPool<RBox_data>& pool)
{
	while (!R.empty())
	{
		pool.put(const_cast<RBox&>(R.top()).move_data());
		R.pop();
	}

	bool ret = true;
	// open file
	std::ifstream fin {fn, std::ios_base::in | std::ios_base::binary};

	if (!fin)
	{
		ret = false;
		goto restore_exit;
	}

	fin.exceptions(std::ifstream::failbit);

	try
	{
		// read basic data
		fin.read(reinterpret_cast<char *>(&optval), sizeof(optval));
		fin.read(reinterpret_cast<char*>(&xopt), sizeof(xopt));
		fin.read(reinterpret_cast<char*>(&runtime), sizeof(runtime));
		fin.read(reinterpret_cast<char*>(&iter), sizeof(iter));
		fin.read(reinterpret_cast<char*>(&lastUpdate), sizeof(lastUpdate));

		// read PBox
		fin.read(reinterpret_cast<char*>(&P), sizeof(P));

		// restore queue
		R.restore(fin, pool);
	}
	catch (...)
	{
		ret = false;
	}

	fin.close();

restore_exit:
	if (!ret)
	{
		std::remove(fn.c_str());

		if (fn.size() < 4 || fn.substr(fn.size()-4,4).compare(".bak") != 0)
			return restore(fn + ".bak", R, P, pool);
	}
	else
		std::cout << ">>> Successfully restored backup" << std::endl;

	return ret;
}

template <size_t Dim>
bool
SIT<Dim>::RType::backup(std::ofstream& fout)
{
	const auto &c = RType_base::c;

	const size_t len = c.size();
	fout.write(reinterpret_cast<const char*>(&len), sizeof(len));

	for (auto it = c.begin(); it != c.end(); ++it)
	{
		fout.write(reinterpret_cast<const char*>(&it->beta), sizeof(it->beta));
		fout.write(reinterpret_cast<const char*>(&it->lb()), sizeof(it->lb()));
		fout.write(reinterpret_cast<const char*>(&it->ub()), sizeof(it->ub()));
		fout.write(reinterpret_cast<const char*>(&it->xk()), sizeof(it->xk()));

		ykptr p = ykptr::undef;

		if (it->yk() == &it->xk())
			p = ykptr::xk;
		else if (it->yk() == &it->ub())
			p = ykptr::ub;
		else if (it->yk() == &it->lb())
			p = ykptr::lb;
		else if (it->yk() == nullptr)
			p = ykptr::null;
		else
		{
			std::cerr << "Invalid yk. Aborting checkpoint." << std::endl;
			return false;
		}

		fout.write(reinterpret_cast<const char*>(&p), sizeof(p));
	}

	return true;
}

template <size_t Dim>
bool
SIT<Dim>::RType::restore(std::ifstream& fin, MiniPool<RBox_data>& pool)
{
	auto &c = RType_base::c;

	// get record count
	size_t len;
	fin.read(reinterpret_cast<char *>(&len), sizeof(len));

	// reserve space
	c.reserve(len);

	for (size_t i = 0; i < len; ++i)
	{
		RBox r(pool.get());

		fin.read(reinterpret_cast<char *>(&r.beta), sizeof(r.beta));
		fin.read(reinterpret_cast<char *>(&r.lb()), sizeof(r.lb()));
		fin.read(reinterpret_cast<char *>(&r.ub()), sizeof(r.ub()));
		fin.read(reinterpret_cast<char *>(&r.xk()), sizeof(r.xk()));

		ykptr p;
		fin.read(reinterpret_cast<char*>(&p), sizeof(p));

		switch (p)
		{
			case ykptr::xk:
				r.yk() = &r.xk();
				break;
				
			case ykptr::ub:
				r.yk() = &r.ub();
				break;
				
			case ykptr::lb:
				r.yk() = &r.lb();
				break;

			case ykptr::null:
				r.yk() = nullptr;
				break;

			default:
				return false;
				break;
		}

		this->push(std::move(r));
	}

	return true;
}


template <size_t Dim>
void
SIT<Dim>::enableBackup(const std::string& filename, const std::chrono::seconds& interval)
{
	backupFile = filename;
	backupInterval = interval;
	backupEnable = true;
}


template <size_t Dim>
void
SIT<Dim>::initBackup()
{
	if (backupEnable)
		backupNext = backupClock::now() + backupInterval;
}


template <size_t Dim>
void
SIT<Dim>::removeBackup() const
{
	std::remove(backupFile.c_str());
	
	const std::string fn = backupFile + ".bak";
	std::remove(fn.c_str());
}

#endif
