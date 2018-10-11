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

#ifndef _MOSEK_H
#define _MOSEK_H

#include <memory>

extern "C"
{
	#include "mosek.h"
	#include "scopt-ext.h"
}

class MSK_error : public std::runtime_error
{
	public:
		MSK_error(const MSKrescodee err, const std::string str);
		const MSKrescodee rescode;
};

void MSKAPI printstr(void *handle, const char str[]);
void mskcall(const MSKrescodee res);

class MosekEnv
{
	public:
		static MSKenv_t getEnv();

		MosekEnv(MosekEnv const&) = delete;
		void operator=(MosekEnv const&) = delete;

		~MosekEnv() noexcept;

	private:
		MosekEnv();
		MSKenv_t env;
		void destroy() noexcept;
};
#endif
