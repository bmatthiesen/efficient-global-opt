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
		typedef std::array<Complex, dimIn> chan_t;
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

		static const size_t wp_drop {1};
		static const size_t wp_snr {0};

		hid_t hComplex_f;
		hid_t hComplex_m ;
		hid_t hChan_f;
		hid_t hWP_m;
		hid_t hWP_f;
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
		void createWPFile(const hsize_t chan_count, const chan_t* h, const chan_t* g, const hsize_t P_size, const double* P, const double NN) const;

		void setResultFile(const std::string outfileBase, const size_t dim1, const size_t dim2 = 0);
		void saveResult(const Result& res, const size_t idx1 = 0, const size_t idx2 = 0);
		void saveResult(Result& res, const size_t idx1 = 0, const size_t idx2 = 0);

		template <size_t NumConstraints>
		void saveResult(const RR<dim,dimC,NumConstraints> * const rr, const size_t idx1 = 0, const size_t idx2 = 0);
		size_t getStartIdx();

	private:
		void createResultFile();
		void emptyDataset(const hid_t loc_id, const char *name, const std::vector<hsize_t> dims, const hid_t hChan_f, const hid_t hChan_m, const void *fill) const;
		void getWPidx(const hsize_t idx, const hid_t file, hsize_t *wpidx) const;
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
	vecIn_t P;
	double P0;
	double N;
	double N0;
	chan_t h;
	chan_t g;
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
	const hsize_t countIn[] = {dimIn};

	// Complex
	hComplex_f = h5call(H5Tcreate(H5T_COMPOUND, sizeof(Complex)));
	h5call(H5Tinsert(hComplex_f, "re", HOFFSET(Complex, re), dtype));
	h5call(H5Tinsert(hComplex_f, "im", HOFFSET(Complex, im), dtype));

	hComplex_m = h5call(H5Tcreate(H5T_COMPOUND, sizeof(Complex)));
	h5call(H5Tinsert(hComplex_m, "re", HOFFSET(Complex, re), H5T_NATIVE_DOUBLE));
	h5call(H5Tinsert(hComplex_m, "im", HOFFSET(Complex, im), H5T_NATIVE_DOUBLE));

	// chan_t
	hChan_f = h5call(H5Tarray_create(hComplex_f, 1, countIn));
	hChan_m = h5call(H5Tarray_create(hComplex_m, 1, countIn));

	// wp indices
	hWP_m = h5call(H5Tcreate(H5T_COMPOUND, 2*sizeof(hsize_t)));
	h5call(H5Tinsert(hWP_m, "channel index", wp_drop * sizeof(hsize_t), H5T_NATIVE_HSIZE));
	h5call(H5Tinsert(hWP_m, "P index", wp_snr * sizeof(hsize_t), H5T_NATIVE_HSIZE));

	hWP_f = h5call(H5Tcreate(H5T_COMPOUND, 2*sizeof(hsize_t)));
	h5call(H5Tinsert(hWP_f, "channel index", 0, H5T_STD_U64LE));
	h5call(H5Tinsert(hWP_f, "P index", sizeof(hsize_t), H5T_STD_U64LE));

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
	H5Tclose(hChan_f);
	H5Tclose(hChan_m);
	H5Tclose(hWP_f);
	H5Tclose(hWP_m);
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

	/*
	 * get WP idx
	 */
	hsize_t wpidx[2];
	getWPidx(WPidx, file, wpidx);

	/*
	 * read WP
	 */
	auto wp = std::make_unique<WP>();

	hid_t g_input = h5call(H5Gopen(file, "/input", H5P_DEFAULT));
	H5Gclose(g_input);
	g_input = h5call(H5Gopen(file, "/input", H5P_DEFAULT));

	// N & N0
	hid_t dataset = h5call(H5Dopen(g_input, "N", H5P_DEFAULT));
	h5call(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(wp->N)));
	H5Dclose(dataset);

	wp->N0 = wp->N;

	// P & P0
	dataset = h5call(H5Dopen(g_input, "P", H5P_DEFAULT));
	hid_t dataspace = h5call(H5Dget_space(dataset));

	const hsize_t count[] = {1};
	h5call(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &wpidx[wp_snr], NULL, count, NULL));

	hsize_t offset[] = {WPidx};
	offset[0] = 1;
	hid_t s_scalar = h5call(H5Screate_simple(1, offset, NULL));

	offset[0] = 0;
	h5call(H5Sselect_hyperslab(s_scalar, H5S_SELECT_SET, offset, NULL, count, NULL));

	h5call(H5Dread(dataset, H5T_NATIVE_DOUBLE, s_scalar, dataspace, H5P_DEFAULT, &(wp->P0)));

	H5Sclose(dataspace);
	H5Dclose(dataset);

	for (hsize_t i = 0; i < dim; ++i)
		wp->P[i] = wp->P0;

	// channels
	dataset = h5call(H5Dopen(g_input, "h", H5P_DEFAULT));
	dataspace = h5call(H5Dget_space(dataset));

	h5call(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &wpidx[wp_drop], NULL, count, NULL));

	h5call(H5Dread(dataset, hChan_m, s_scalar, dataspace, H5P_DEFAULT, &(wp->h)));
	H5Dclose(dataset);
	H5Sclose(dataspace);

	dataset = h5call(H5Dopen(g_input, "g", H5P_DEFAULT));
	dataspace = h5call(H5Dget_space(dataset));

	h5call(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &wpidx[wp_drop], NULL, count, NULL));

	h5call(H5Dread(dataset, hChan_m, s_scalar, dataspace, H5P_DEFAULT, &(wp->g)));
	H5Dclose(dataset);
	H5Sclose(dataspace);

	H5Gclose(g_input);
	H5Sclose(s_scalar);
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


template <size_t dim, size_t dimC, size_t dimIn>
void
HDF<dim, dimC, dimIn>::createWPFile(const hsize_t chan_count, const chan_t* h, const chan_t* g, const hsize_t P_size, const double* P, const double N) const
{
	hid_t file = h5call(H5Fcreate(wpfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)); // TODO H5F_ACC_EXCL

	/*
	 * input data
	 */
	hid_t g_input = h5call(H5Gcreate(file, "/input", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5Gclose(g_input);
	g_input = h5call(H5Gopen(file, "/input", H5P_DEFAULT));

	hid_t dataspace, dataset;

	// N
	hsize_t count[] = {1};
	dataspace = h5call(H5Screate_simple(1, count, NULL));

	dataset = h5call(H5Dcreate(g_input, "N", dtype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	h5call(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N));
	H5Dclose(dataset);

	// P
	h5call(H5Sset_extent_simple(dataspace, 1, &P_size, NULL));

	dataset = h5call(H5Dcreate(g_input, "P", dtype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	h5call(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, P));
	H5Dclose(dataset);

	// channels
	h5call(H5Sset_extent_simple(dataspace, 1, &chan_count, NULL));

	dataset = h5call(H5Dcreate(g_input, "h", hChan_f, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	h5call(H5Dwrite(dataset, hChan_m, H5S_ALL, H5S_ALL, H5P_DEFAULT, h));
	H5Dclose(dataset);

	dataset = h5call(H5Dcreate(g_input, "g", hChan_f, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	h5call(H5Dwrite(dataset, hChan_m, H5S_ALL, H5S_ALL, H5P_DEFAULT, g));
	H5Dclose(dataset);

	// close group
	h5call(H5Gclose(g_input));

	/*
	 * WPs
	 */
	const hsize_t numWP = chan_count * P_size;
	hsize_t wp[numWP][2];

	size_t idx = 0;
	for (hsize_t drop = 0; drop < chan_count; ++drop)
	{
		for (hsize_t snr = 0; snr < P_size; ++snr)
		{
			wp[idx][wp_drop] = drop;
			wp[idx][wp_snr] = snr;
			++idx;
		}
	}

	h5call(H5Sset_extent_simple(dataspace, 1, &numWP, NULL));

	dataset = h5call(H5Dcreate(file, "WP", hWP_f, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	h5call(H5Dwrite(dataset, hWP_m, H5S_ALL, H5S_ALL, H5P_DEFAULT, wp));
	H5Dclose(dataset);

	/*
	 * Clean up
	 */
	H5Sclose(dataspace);
	H5Fclose(file);
}


template <size_t dim, size_t dimC, size_t dimIn>
void
HDF<dim, dimC, dimIn>::getWPidx(const hsize_t idx, const hid_t file, hsize_t *wpidx) const
{

	hid_t dataset = h5call(H5Dopen(file, "WP", H5P_DEFAULT));
	hid_t dataspace = h5call(H5Dget_space(dataset));

	const hsize_t count[] = {1};
	hsize_t offset[] = {idx};
	h5call(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL));

	offset[0] = 1;
	hid_t s_scalar = h5call(H5Screate_simple(1, offset, NULL));

	offset[0] = 0;
	h5call(H5Sselect_hyperslab(s_scalar, H5S_SELECT_SET, offset, NULL, count, NULL));

	h5call(H5Dread(dataset, hWP_m, s_scalar, dataspace, H5P_DEFAULT, wpidx));

	H5Sclose(dataspace);
	H5Dclose(dataset);
}
#endif
