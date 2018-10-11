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

#include <string>
#include <stdexcept>

extern "C" {
	#include "hdf5.h"
}

#include "HDF.h"

hid_t
h5call(const hid_t ret)
{
	if (ret < 0)
		throw std::runtime_error("H5 returned error: " + std::to_string(ret)); // error

	return ret;
}

bool
fexist (const std::string& name) {
	struct stat buffer;   
	return (stat (name.c_str(), &buffer) == 0); 
}


std::string
create_filename(const std::string& fname, const std::string& prefix, const std::string& suffix)
{
	// find basename
	const std::string::size_type end = fname.rfind(".");
	std::string::size_type start = fname.rfind("/");

	std::string ret;

	if (start == std::string::npos)
		start = 0;
	else
	{
		++start;
		ret += fname.substr(0,start);
	}

	ret += prefix;

	if (end < fname.size()-1 && end >= start)
		ret += fname.substr(start, end-start) + suffix + fname.substr(end);
	else
		ret += fname.substr(start) + suffix;

	return ret;
}
