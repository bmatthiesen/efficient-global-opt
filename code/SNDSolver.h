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

#ifndef _SNDSOLVER_H
#define _SNDSOLVER_H

#include <memory>
#include <iostream>
#include <sstream>
using std::cout; using std::endl;

#include "SND.h"
#include "HDF.h"

template <size_t C = 8, bool useTighterBound = true>
class SNDSolver
{
	static_assert(C <= 7, "Out of range");

public:
	SND<C> s;

	using H5 = HDF<SNDDIM,SNDDIM>;

	template <typename... Ts>
	SNDSolver(H5::WPPtr wp, Ts&&... args) try
		: s(std::forward<Ts>(args)...)
	{
		s.init(wp, useTighterBound);

		typename SND<C>::vtypeC linC;
		linC.fill(1.0);
		
		s.setObjective(linC, {});
	}
	catch (GRBException &e)
	{
		std::stringstream ss;
		ss << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")";
		throw std::runtime_error(ss.str());
	}

	virtual ~SNDSolver() {};

	void setObjective(typename SND<0>::vtypeC linC, typename SND<0>::vtype linNC, typename SND<0>::vtypeC denomLinC = {}, typename SND<0>::vtype denomLinNC = {}, const basetype numerConst = 0.0, const basetype denomConst = 0.0)
	{
		s.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
	}

	virtual std::string getName() const
	{
		if constexpr(C == 7)
			return "afsnd";
		else if constexpr(C == 0)
			return "aftin";
		else
			return std::string("af") + std::to_string(C);
	}

	constexpr size_t getDim()
	{
		return 1;
	}

	const SND<C>* solve(const size_t idx)
	{
		if (idx >= getDim())
			throw std::runtime_error("foo");

		try
		{
			s.optimize();
		}
		catch (GRBException &e)
		{
			std::stringstream ss;
			ss << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")";
			throw std::runtime_error(ss.str());
		}

		return &s;
	}

	void setPrecision(const basetype eta, const basetype epsilon)
	{
		s.setPrecision(eta, epsilon);
	}

	void enableBackup(const std::string& filename, const std::chrono::seconds& interval)
	{
		s.enableBackup(filename, interval);
	}

	void removeBackup() const
	{
		s.removeBackup();
	}
};

template <bool useTighterBound>
class SNDSolver<8, useTighterBound>
{
public:
	SND<0> s0;
	SND<1> s1;
	SND<2> s2;
	SND<3> s3;
	SND<4> s4;
	SND<5> s5;
	SND<6> s6;
	SND<7> s7;

	using H5 = HDF<SNDDIM,SNDDIM>;

	template <typename... Ts>
	SNDSolver(H5::WPPtr wp, Ts&&... args)
		: 	s0(std::forward<Ts>(args)...),
			s1(std::forward<Ts>(args)...),
			s2(std::forward<Ts>(args)...),
			s3(std::forward<Ts>(args)...),
			s4(std::forward<Ts>(args)...),
			s5(std::forward<Ts>(args)...),
			s6(std::forward<Ts>(args)...),
			s7(std::forward<Ts>(args)...)
	{
		s0.init(wp, useTighterBound);
		s1.init(wp, useTighterBound);
		s2.init(wp, useTighterBound);
		s3.init(wp, useTighterBound);
		s4.init(wp, useTighterBound);
		s5.init(wp, useTighterBound);
		s6.init(wp, useTighterBound);
		s7.init(wp, useTighterBound);

		typename SND<0>::vtypeC linC;
		linC.fill(1.0);
		setObjective(linC, {});
	}

	virtual ~SNDSolver() {};

	void setObjective(typename SND<0>::vtypeC linC, typename SND<0>::vtype linNC, typename SND<0>::vtypeC denomLinC = {}, typename SND<0>::vtype denomLinNC = {}, const basetype numerConst = 0.0, const basetype denomConst = 0.0)
	{
		s0.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
		s1.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
		s2.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
		s3.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
		s4.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
		s5.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
		s6.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
		s7.setObjective(linC, linNC, denomLinC, denomLinNC, numerConst, denomConst);
	}

	void setPrecision(const basetype eta, const basetype epsilon)
	{
		s0.setPrecision(eta, epsilon);
		s1.setPrecision(eta, epsilon);
		s2.setPrecision(eta, epsilon);
		s3.setPrecision(eta, epsilon);
		s4.setPrecision(eta, epsilon);
		s5.setPrecision(eta, epsilon);
		s6.setPrecision(eta, epsilon);
		s7.setPrecision(eta, epsilon);
	}

	virtual std::string getName() const
	{
		return "afsnd2";
	}

	constexpr size_t getDim()
	{
		return 8;
	}

	const H5::Result solve(const size_t idx)
	{
		try
		{
			switch (idx) {
			case 0:
				s0.optimize();
				{
					H5::Result res(&s0);
					return res;
				}
			case 1:
				s1.optimize();
				{
					H5::Result res(&s1);
					return res;
				}
			case 2:
				s2.optimize();
				{
					H5::Result res(&s2);
					return res;
				}
			case 3:
				s3.optimize();
				{
					H5::Result res(&s3);
					return res;
				}
			case 4:
				s4.optimize();
				{
					H5::Result res(&s4);
					return res;
				}
			case 5:
				s5.optimize();
				{
					H5::Result res(&s5);
					return res;
				}
			case 6:
				s6.optimize();
				{
					H5::Result res(&s6);
					return res;
				}
			case 7:
				s7.optimize();
				{
					H5::Result res(&s7);
					return res;
				}
			default:
				throw std::runtime_error("idx too big");
			}
		}
		catch (GRBException &e)
		{
			std::stringstream ss;
			ss << "Gurobi exception: " << e.getMessage() << "(" << e.getErrorCode() << ")";
			throw std::runtime_error(ss.str());
		}
	}

	void enableBackup(const std::string& filename, const std::chrono::seconds& interval)
	{
		s0.enableBackup(filename + "0", interval);
		s1.enableBackup(filename + "1", interval);
		s2.enableBackup(filename + "2", interval);
		s3.enableBackup(filename + "3", interval);
		s4.enableBackup(filename + "4", interval);
		s5.enableBackup(filename + "5", interval);
		s6.enableBackup(filename + "6", interval);
		s7.enableBackup(filename + "7", interval);
	}

	void removeBackup() const
	{
		s0.removeBackup();
		s1.removeBackup();
		s2.removeBackup();
		s3.removeBackup();
		s4.removeBackup();
		s5.removeBackup();
		s6.removeBackup();
		s7.removeBackup();
	}
};

#endif
