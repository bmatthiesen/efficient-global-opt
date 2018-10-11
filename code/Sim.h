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

#ifndef _SIM_H
#define _SIM_H

#include <string>
#include <sstream>
#include <cstdlib>
#include <chrono>

#include "gurobi_c++.h"

#include <type_traits>
#include <memory>
template <class T> struct is_smart_ptr : std::false_type {};
template <class T> struct is_smart_ptr<std::unique_ptr<T>> : std::true_type {};
template <class T> struct is_smart_ptr<std::shared_ptr<T>> : std::true_type {};

constexpr std::chrono::hours backupInterval(5);

template <class S>
class Sim
{
	protected:
		typename S::H5 hdf;
		typename S::H5::WPPtr wp;
		S solver;


	private:
#if defined(SLURM) && defined(SLURMLOCK)
		std::string lockfile;
#endif

		struct argparse;
		
		template <typename... Ts>
		Sim(const argparse& args, Ts&&... fwd);

#if defined(SLURM) && defined(SLURMLOCK)
		void lock() const;
		void unlock() const;
#endif


	public:
		template <typename... Ts>
		Sim(int argc, char *argv[], Ts&&... args)
			: Sim(argparse(argc, argv), std::forward<Ts>(args)...)
		{}

		virtual ~Sim() {}

		virtual void run();
};


template <class S>
struct Sim<S>::argparse
{
	size_t WPidx;
	std::string wpfile;

	argparse(int argc, char *argv[]);
};



template <class S>
void
Sim<S>::run()
{
#if defined(SLURM) && defined(SLURMLOCK)
	lock();
#endif

	size_t start = hdf.getStartIdx();

	if (start > 0)
	{
		std::cout << "Result file exists: ";

		if (start < solver.getDim())
			std::cout << "Starting with index " << start << std::endl;
		else
			std::cout << "Nothing left to do" << std::endl;
	}

	for (size_t idx = start; idx < solver.getDim(); ++idx)
	{
		auto r = solver.solve(idx);

		if constexpr(is_smart_ptr<decltype(r)>::value)
			hdf.saveResult(*r.get(), idx);
		else
			hdf.saveResult(r, idx);
	}

	solver.removeBackup();

#if defined(SLURM) && defined(SLURMLOCK)
	unlock();
#endif
}




template <class S>
template <typename... Ts>
Sim<S>::Sim(const argparse& args, Ts&&... fwd)
	: hdf(args.wpfile, args.WPidx), wp(hdf.getWP()), solver(wp, std::forward<Ts>(fwd)...)
{
#ifdef SLURM
	const char * const e = getenv("JOB_HPC_SAVEDIR");

	if (e == nullptr)
		throw std::runtime_error("can't get savedir from JOB_HPC_SAVEDIR");

#ifdef SLURMLOCK
	lockfile = std::string(e) + "/" + solver.getName() + "_wp" + std::to_string(args.WPidx) + ".lock";
#endif

	hdf.setResultFile(std::string(e) + "/res_" + solver.getName() + "_", solver.getDim());

	solver.enableBackup(std::string(e) + "/ckpt_" + solver.getName() + "_wp=" + std::to_string(args.WPidx), backupInterval);
#else
	hdf.setResultFile("../results/res_" + solver.getName() + "_", solver.getDim());
	solver.enableBackup("../results/ckpt_" + solver.getName() + "_wp=" + std::to_string(args.WPidx), backupInterval);
#endif
}

template <class S>
Sim<S>::argparse::argparse(int argc, char *argv[])
{
#ifdef SLURM
	if (argc != 2)
#else
	if (argc != 3)
#endif
	{
		std::stringstream ss;
		ss << argv[0] << ": Wrong number of arguments: " << argv[0] << " wpfile";

#ifndef SLURM
		ss << " wpindex";
#endif

		throw std::runtime_error(ss.str());
	}

	wpfile = argv[1];

#ifdef SLURM
	char *e = getenv("SLURM_ARRAY_TASK_ID");

	if (e == nullptr)
		throw std::runtime_error("can't get WP idx from SLURM_ARRAY_TASK_ID");

	WPidx = strtoull(e, NULL, 10);
#else
	WPidx = strtoull(argv[2], NULL, 10);
#endif
}


#if defined(SLURM) && defined(SLURMLOCK)
template <class S>
void
Sim<S>::lock() const
{
	FILE *fp;

	fp = std::fopen(lockfile.c_str(),"r");

	if (fp)
	{
		std::fclose(fp);
		throw std::runtime_error("lock present");
	}

	fp = std::fopen(lockfile.c_str(), "w");
	std::fclose(fp);
}

template <class S>
inline void
Sim<S>::unlock() const
{
	std::remove(lockfile.c_str());
}
#endif

#endif
