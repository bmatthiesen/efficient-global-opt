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

#ifndef _HDF_H
#define _HDF_H

#include <memory>
#include <string>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>

#include <sys/stat.h>

extern "C" {
	#include <unistd.h>
	#include "hdf5.h"
}

#include "SIT.h"
#include "RR.h"
#include "util.h"

using std::size_t;

hid_t h5call(const hid_t ret);
bool fexist (const std::string& name);
std::string create_filename(const std::string& path, const std::string& prefix, const std::string& suffix);

template <size_t dim, size_t dimC, size_t dimIn = dim>
class HDF
{
	public:
		struct WP;
		struct CP;
		struct Result;
		struct Complex;
		typedef std::shared_ptr<const WP> WPPtr;
		typedef Complex chan_t[dimIn][dimIn];
		typedef std::array<double, dimIn> vecIn_t;
		typedef std::array<double, dim> vec_t;
		typedef std::array<double, dimC> vecC_t;

		enum class SolverStatus { Optimal, Unsolved, Infeasible, Maxiter };

		hid_t hChan_m;

		const hsize_t WPidx;
	private:
		const std::string wpfile;
		std::string outfile;
		const hid_t dtype;

		size_t dim1, dim2;
		bool RFexists;

		hid_t hComplex_f;
		hid_t hComplex_m ;
		hid_t hPopt_f;
		hid_t hPopt_m;
		hid_t hxoptC_f;
		hid_t hxoptC_m;
		hid_t hStatus;
		hid_t hResult_m;
		hid_t hResult_f;
		hid_t hCP_m;
		hid_t hCP_f;

	public:
		HDF(const std::string wpfile, const hsize_t WPidx, const hid_t double_type = H5T_IEEE_F64LE);
		~HDF();
		std::unique_ptr<WP> getWP() const;

		void setResultFile(const std::string outfileBase, const size_t dim1, const size_t dim2 = 0);
		void saveResult(const Result& res, const size_t idx1 = 0, const size_t idx2 = 0);
		void saveResult(Result& res, const size_t idx1 = 0, const size_t idx2 = 0);

		template <size_t NumConstraints>
		void saveResult(const RR<dim,dimC,NumConstraints> * const rr, const size_t idx1 = 0, const size_t idx2 = 0);
		size_t getStartIdx();

	private:
		void createResultFile();
};

template <size_t dim, size_t dimC, size_t dimIn>
struct HDF<dim, dimC, dimIn>::Complex
{
	double re;
	double im;
};

template <size_t dim, size_t dimC, size_t dimIn>
struct HDF<dim, dimC, dimIn>::WP
{
	chan_t h;
};

template <size_t dim, size_t dimC, size_t dimIn>
struct HDF<dim, dimC, dimIn>::Result
{
	hsize_t WPidx;
	double objective;
	vecC_t xopt_c;
	vec_t xopt_nc;
	SolverStatus status;
	double runtime;
	unsigned long long iter;
	unsigned long long lastUpdate;
	hsize_t peakRSS;

	Result();

	template <size_t NumConstraints>
	Result(const RR<dim,dimC,NumConstraints> * const rr);
};

template <size_t dim, size_t dimC, size_t dimIn>
inline
HDF<dim, dimC, dimIn>::Result::Result()
	: WPidx(0), objective(0), xopt_c {}, xopt_nc {}, status(SolverStatus::Unsolved), runtime(0), iter(0), lastUpdate(0), peakRSS(getPeakRSS())
{ }

template <size_t dim, size_t dimC, size_t dimIn>
template <size_t NumConstraints>
HDF<dim, dimC, dimIn>::Result::Result(const RR<dim,dimC,NumConstraints> * const rr)
	: WPidx(0), objective(rr->optval), runtime(rr->runtime), iter(rr->iter), lastUpdate(rr->lastUpdate), peakRSS(getPeakRSS())
{
	std::copy_n(rr->optNC.begin(), xopt_nc.size(), xopt_nc.begin());
	std::copy_n(rr->optC.begin(), xopt_c.size(), xopt_c.begin());

	switch (rr->status) {
	case RR<dim,dimC,NumConstraints>::Status::Optimal:
		status = SolverStatus::Optimal;
		break;

	case RR<dim,dimC,NumConstraints>::Status::Infeasible:
		status = SolverStatus::Infeasible;
		break;

	case RR<dim,dimC,NumConstraints>::Status::Unsolved:
		status = SolverStatus::Unsolved;
		break;
	}

}

template <size_t dim, size_t dimC, size_t dimIn>
struct HDF<dim, dimC, dimIn>::CP
{
	vec_t Popt;
	double objective;
	unsigned long long iter;
	double runtime;
};


template <size_t dim, size_t dimC, size_t dimIn>
HDF<dim, dimC, dimIn>::HDF(const std::string wpfile, const hsize_t WPidx, const hid_t double_type)
	: WPidx(WPidx), wpfile(wpfile), dtype(double_type), dim1(0), dim2(0), RFexists(false)
{
	/* Datatypes */
	const hsize_t count[] = {dim};
	const hsize_t countC[] = {dimC};

	// Complex
	hComplex_f = h5call(H5Tcreate(H5T_COMPOUND, sizeof(Complex)));
	h5call(H5Tinsert(hComplex_f, "r", HOFFSET(Complex, re), dtype));
	h5call(H5Tinsert(hComplex_f, "i", HOFFSET(Complex, im), dtype));

	hComplex_m = h5call(H5Tcreate(H5T_COMPOUND, sizeof(Complex)));
	h5call(H5Tinsert(hComplex_m, "r", HOFFSET(Complex, re), H5T_NATIVE_DOUBLE));
	h5call(H5Tinsert(hComplex_m, "i", HOFFSET(Complex, im), H5T_NATIVE_DOUBLE));

	// Popt
	hPopt_f = h5call(H5Tarray_create(dtype, 1, count));
	hPopt_m = h5call(H5Tarray_create(H5T_NATIVE_DOUBLE, 1, count));

	hxoptC_f = h5call(H5Tarray_create(dtype, 1, countC));
	hxoptC_m = h5call(H5Tarray_create(H5T_NATIVE_DOUBLE, 1, countC));

	// Status
	hStatus = h5call(H5Tenum_create(H5T_STD_U16BE));

	std::uint16_t val = static_cast<std::uint16_t>(SolverStatus::Optimal);
	h5call(H5Tenum_insert(hStatus, "Optimal", &val));

	val = static_cast<std::uint16_t>(SolverStatus::Infeasible);
	h5call(H5Tenum_insert(hStatus, "Infeasible", &val));

	val = static_cast<std::uint16_t>(SolverStatus::Maxiter);
	h5call(H5Tenum_insert(hStatus, "Maxiter", &val));

	val = static_cast<std::uint16_t>(SolverStatus::Unsolved);
	h5call(H5Tenum_insert(hStatus, "Unsolved", &val));

	// Result
	hResult_m = h5call(H5Tcreate(H5T_COMPOUND, sizeof(Result)));
	h5call(H5Tinsert(hResult_m, "WP index", HOFFSET(Result, WPidx), H5T_NATIVE_HSIZE));
	h5call(H5Tinsert(hResult_m, "Objective Value", HOFFSET(Result, objective), H5T_NATIVE_DOUBLE));
	h5call(H5Tinsert(hResult_m, "xopt_C", HOFFSET(Result, xopt_c), hxoptC_m));
	h5call(H5Tinsert(hResult_m, "xopt_NC", HOFFSET(Result, xopt_nc), hPopt_m));
	h5call(H5Tinsert(hResult_m, "Status", HOFFSET(Result, status), hStatus));
	h5call(H5Tinsert(hResult_m, "Runtime", HOFFSET(Result, runtime), H5T_NATIVE_DOUBLE));
	h5call(H5Tinsert(hResult_m, "Iterations", HOFFSET(Result, iter), H5T_NATIVE_ULLONG));
	h5call(H5Tinsert(hResult_m, "Last Update", HOFFSET(Result, lastUpdate), H5T_NATIVE_ULLONG));
	h5call(H5Tinsert(hResult_m, "Peak RSS", HOFFSET(Result, peakRSS), H5T_NATIVE_HSIZE));

	hResult_f = h5call(H5Tcreate(H5T_COMPOUND, sizeof(Result)));
	h5call(H5Tinsert(hResult_f, "WP index", HOFFSET(Result, WPidx), H5T_STD_U64LE));
	h5call(H5Tinsert(hResult_f, "Objective Value", HOFFSET(Result, objective), dtype));
	h5call(H5Tinsert(hResult_f, "xopt_C", HOFFSET(Result, xopt_c), hxoptC_f));
	h5call(H5Tinsert(hResult_f, "xopt_NC", HOFFSET(Result, xopt_nc), hPopt_f));
	h5call(H5Tinsert(hResult_f, "Status", HOFFSET(Result, status), hStatus));
	h5call(H5Tinsert(hResult_f, "Runtime", HOFFSET(Result, runtime), dtype));
	h5call(H5Tinsert(hResult_f, "Iterations", HOFFSET(Result, iter), H5T_STD_U64LE));
	h5call(H5Tinsert(hResult_f, "Last Update", HOFFSET(Result, lastUpdate), H5T_STD_U64LE));
	h5call(H5Tinsert(hResult_f, "Peak RSS", HOFFSET(Result, peakRSS), H5T_STD_U64LE));

	// CP
	hCP_m = h5call(H5Tcreate(H5T_COMPOUND, sizeof(CP)));
	h5call(H5Tinsert(hCP_m, "P_opt", HOFFSET(CP, Popt), hPopt_m));
	h5call(H5Tinsert(hCP_m, "Objective", HOFFSET(CP, objective), H5T_NATIVE_DOUBLE));
	h5call(H5Tinsert(hCP_m, "Iterations", HOFFSET(CP, iter), H5T_NATIVE_ULLONG));
	h5call(H5Tinsert(hCP_m, "Runtime", HOFFSET(CP, runtime), H5T_NATIVE_DOUBLE));

	hCP_f = h5call(H5Tcreate(H5T_COMPOUND, sizeof(CP)));
	h5call(H5Tinsert(hCP_f, "P_opt", HOFFSET(CP, Popt), hPopt_m));
	h5call(H5Tinsert(hCP_f, "Objective", HOFFSET(CP, objective), dtype));
	h5call(H5Tinsert(hCP_f, "Iterations", HOFFSET(CP, iter), H5T_STD_U64LE));
	h5call(H5Tinsert(hCP_f, "Runtime", HOFFSET(CP, runtime), dtype));
}

template <size_t dim, size_t dimC, size_t dimIn>
HDF<dim, dimC, dimIn>::~HDF()
{
	H5Tclose(hComplex_f);
	H5Tclose(hComplex_m);
	H5Tclose(hPopt_f);
	H5Tclose(hPopt_m);
	H5Tclose(hxoptC_f);
	H5Tclose(hxoptC_m);
	H5Tclose(hStatus);
	H5Tclose(hResult_m);
	H5Tclose(hResult_f);
	H5Tclose(hCP_m);
	H5Tclose(hCP_f);
}

template <size_t dim, size_t dimC, size_t dimIn>
std::unique_ptr<typename HDF<dim, dimC, dimIn>::WP>
HDF<dim, dimC, dimIn>::getWP() const
{
	std::cout << "Reading WP " << WPidx << " from file '" << wpfile << "'" << std::endl;

	hid_t file = h5call(H5Fopen(wpfile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
	auto wp = std::make_unique<WP>();

	const hsize_t count[] = {1, dimIn, dimIn};
	hid_t s_channel = h5call(H5Screate_simple(3, count, NULL));

	// channels
	hid_t dataset = h5call(H5Dopen(file, "channel", H5P_DEFAULT));
	hid_t dataspace = h5call(H5Dget_space(dataset));

	const hsize_t index[] = {WPidx, 0, 0};
	h5call(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, index, NULL, count, NULL));

	h5call(H5Dread(dataset, hComplex_m, s_channel, dataspace, H5P_DEFAULT, &(wp->h)));
	H5Dclose(dataset);
	H5Sclose(dataspace);

	H5Sclose(s_channel);
	H5Fclose(file);

	return wp;
}


template <size_t dim, size_t dimC, size_t dimIn>
void
HDF<dim, dimC, dimIn>::setResultFile(const std::string outfileBase, const size_t dim1, const size_t dim2)
{
	if (RFexists)
		throw std::runtime_error("Result file already created!");

	if (dim1 <= 0)
		throw std::runtime_error("dim1 too small");

	if (dim2 < 0)
		throw std::runtime_error("dim2 too small");

	this->dim1 = dim1;
	this->dim2 = dim2;

	outfile = outfileBase + "wp=" + std::to_string(WPidx) + ".h5";
}


template <size_t dim, size_t dimC, size_t dimIn>
void
HDF<dim, dimC, dimIn>::createResultFile()
{
	if (RFexists)
		throw std::runtime_error("already created result file");

	if (dim1 <= 0)
		throw std::runtime_error("dim1 too small");

	if (dim2 < 0)
		throw std::runtime_error("dim2 too small");

	/*
	 * HDF
	 */
	hid_t file = h5call(H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));

	const hsize_t count[] = {dim1, dim2};
	const hsize_t rank = dim2 == 0 ? 1 : 2;
	hid_t dataspace = h5call(H5Screate_simple(rank, count, NULL));

	hid_t dataset = h5call(H5Dcreate(file, "results", hResult_f, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

	H5Dclose(dataset);
	H5Sclose(dataspace);
	H5Fclose(file);

	RFexists = true;
}

template <size_t dim, size_t dimC, size_t dimIn>
size_t
HDF<dim, dimC, dimIn>::getStartIdx()
{
	if (RFexists)
		throw std::runtime_error("already created result file");

	if (access(outfile.c_str(), R_OK|W_OK) == 0 && H5Fis_hdf5(outfile.c_str()))
	{
		hid_t file = h5call(H5Fopen(outfile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));

		// open dataset
		hid_t dataset = h5call(H5Dopen(file, "results", H5P_DEFAULT));
		hid_t filespace = h5call(H5Dget_space(dataset));

		// sanity chceks
		bool badfile = false;
		if (H5Sget_simple_extent_ndims(filespace) != 1)
			badfile = true;

		hsize_t dim11;
		h5call(H5Sget_simple_extent_dims(filespace, &dim11, nullptr));
		if (dim11 != dim1)
			badfile = true;

		if (badfile)
		{
			H5Dclose(dataset);
			H5Sclose(filespace);
			H5Fclose(file);

			return 0;
		}

		// looks good. keep it
		RFexists = true;

		// singleton memspace
		const hsize_t count[] = {1};
		const hsize_t zero[] = {0};
		hid_t memspace = h5call(H5Screate_simple(1, count, NULL));
		h5call(H5Sselect_hyperslab(memspace, H5S_SELECT_SET, zero, NULL, count, NULL));

		// search for first empty result
		hsize_t i;
		for (i = 0; i < dim1; ++i)
		{
			const hsize_t idx[] = {i};
			Result data;

			h5call(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, idx, NULL, count, NULL));
			h5call(H5Dread(dataset, hResult_m, memspace, filespace, H5P_DEFAULT, &data));

			if (data.WPidx == 0 && data.objective == 0)
				break;
		}

		H5Dclose(dataset);
		H5Sclose(filespace);
		H5Sclose(memspace);
		H5Fclose(file);

		return i;
	}
	else
		return 0;
}

template <size_t dim, size_t dimC, size_t dimIn>
void
HDF<dim, dimC, dimIn>::saveResult(const Result& res, const size_t idx1, const size_t idx2)
{
	if (!RFexists)
		createResultFile();

	if (dim2 == 0)
	{
		if (dim1 == 1)
			std::cout << "Storing result to file '" << outfile << "'" << std::endl;
		else
			std::cout << "Storing result with idx " << idx1 << " to file '" << outfile << "'" << std::endl;
	}
	else
		std::cout << "Storing result with idx (" << idx1 << "," << idx2 << ") to file '" << outfile << "'" << std::endl;

	/*
	 * Error checks
	 */
	if (idx1 >= dim1 || (dim2 > 0 && idx2 >= dim2))
		throw std::runtime_error("idx too large");

	/*
	 * offsets
	 */
	const hsize_t zero[] = {0};
	const hsize_t count[] = {1, 1};
	const hsize_t idx[] = {idx1, idx2};

	/*
	 * HDF
	 */
	hid_t file = h5call(H5Fopen(outfile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));

	// open dataset and select
	hid_t dataset = h5call(H5Dopen(file, "results", H5P_DEFAULT));
	hid_t filespace = h5call(H5Dget_space(dataset));
	h5call(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, idx, NULL, count, NULL));

	// create singleton memspace
	hid_t memspace = h5call(H5Screate_simple(1, count, NULL));
	h5call(H5Sselect_hyperslab(memspace, H5S_SELECT_SET, zero, NULL, count, NULL));
	
	// write
	h5call(H5Dwrite(dataset, hResult_m, memspace, filespace, H5P_DEFAULT, &res));

	// close
	H5Dclose(dataset);
	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Fclose(file);
}


template <size_t dim, size_t dimC, size_t dimIn>
void
HDF<dim, dimC, dimIn>::saveResult(Result& res, const size_t idx1, const size_t idx2)
{
	if (res.WPidx == 0)
		res.WPidx = WPidx;

	saveResult(static_cast<const Result&>(res), idx1, idx2);
}



template <size_t dim, size_t dimC, size_t dimIn>
template <size_t NumConstraints>
void
HDF<dim, dimC, dimIn>::saveResult(const RR<dim,dimC,NumConstraints> * const rr, const size_t idx1, const size_t idx2)
{
	Result res(rr);
	res.WPidx = WPidx;

	saveResult(res, idx1, idx2);
}

#endif
