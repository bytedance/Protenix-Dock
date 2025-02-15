/*
* Copyright (C) 2025 ByteDance and/or its affiliates
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <boost/python.hpp>

#include "bytedock/version.h"
#include "bytedock/core/entrance.h"

namespace bytedock {

char const* version() {
    return VERSION_STR;
}

}

using namespace boost::python;
using namespace bytedock;

BOOST_PYTHON_MODULE(_protenix_dock) {
    def("version", &version);
    def("log2file", &setup_global_logging, (arg("log_file") = "-", arg("verbosity") = 1));
    class_<ReusableEngine>("ReusableEngine", init<const std::string&, const std::string&, int>(
        (arg("receptor_file"), arg("sf_file") = "", arg("cpu_nthreads") = 0)
    ))
        .def("set_box", &ReusableEngine::set_box, (
            arg("center_x"), arg("center_y"), arg("center_z"),
            arg("size_x"), arg("size_y"), arg("size_z")
        ))
        .def("generate_maps", &ReusableEngine::generate_maps, (
            arg("grid_spacing"), args("max_size") = 0.
        ))
        .def("dump_maps", &ReusableEngine::dump_maps, arg("cache_dir"))
        .def("load_maps", &ReusableEngine::load_maps, arg("cache_dir"))
        .def("drop_maps", &ReusableEngine::drop_maps)
        .def("evaluate", &ReusableEngine::evaluate, (
            arg("ligand_index_file"), arg("output_dir")
        ))
        .def("optimize", &ReusableEngine::optimize, (
            arg("ligand_index_file"), arg("output_dir"),
            arg("max_niters") = 240, arg("slope") = 1e6
        ))
        .def("search", &ReusableEngine::search, (
            arg("ligand_index_file"), arg("output_dir"),
            arg("seed") = 0, arg("exhaustiveness") = 256, arg("max_nsteps") = 200,
            arg("relax_nsteps") = 40, arg("num_modes") = 8, arg("min_rmsd") = 0.5,
            arg("mc_mmenergy_threshold") = -1., arg("slope") = 1e6
        ))
    ;
}
