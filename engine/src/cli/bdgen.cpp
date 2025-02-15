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

#include <boost/program_options.hpp>

#include "bytedock/core/entrance.h"
#include "bytedock/ext/pfile.h"
#include "bytedock/lib/utility.h"
#include "bytedock/version.h"

namespace bd = bytedock;
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    po::options_description in_args("Inputs");
    std::string sf_file, receptor_file;
    int nthreads;
    in_args.add_options()
        ("scoring_function", po::value<std::string>(&sf_file), "A Json file defining the scoring function. If empty, v6 is enabled.")
        ("receptor_file", po::value<std::string>(&receptor_file), "A Json file describing the receptor.")
        ("nthreads", po::value<int>(&nthreads)->default_value(0), "Thread count for generating maps. 0 means all available cores.")
    ;

    po::options_description out_args("Outputs");
    std::string output_dir, log_file;
    int verbosity;
    out_args.add_options()
        ("output_dir", po::value<std::string>(&output_dir)->default_value("docked"), "Output folder for maps of each atom type.")
        ("log_file", po::value<std::string>(&log_file)->default_value("-"), "Text file for log messages. \"-\" is short for STDOUT.")
        ("verbosity", po::value<int>(&verbosity)->default_value(1), "Verbosity level. (0=warning, 1=info, 2=debug)")
    ;

    po::options_description search_args("Precalculation paramters");
    double center_x, center_y, center_z;
    double size_x, size_y, size_z, max_size;
    double grid_spacing;
    search_args.add_options()
        ("center_x", po::value<double>(&center_x), "X coordinate of box center (Angstrom).")
        ("center_y", po::value<double>(&center_y), "Y coordinate of box center (Angstrom).")
        ("center_z", po::value<double>(&center_z), "Z coordinate of box center (Angstrom).")
        ("size_x", po::value<double>(&size_x), "Box size in X dimension (Angstrom).")
        ("size_y", po::value<double>(&size_y), "Box size in Y dimension (Angstrom).")
        ("size_z", po::value<double>(&size_z), "Box size in Z dimension (Angstrom).")
        ("max_size", po::value<double>(&max_size)->default_value(0), "max size(=0) to generate, map_size_x/y/z = max(max_size, size_x)")
        ("spacing", po::value<double>(&grid_spacing)->default_value(0.375), "Grid spacing with value>=0.01. (Angstrom)")
    ;

    po::options_description config_args("Configuration file (optional)");
    std::string config_file;
    config_args.add_options()
        ("config", po::value<std::string>(&config_file), "The above options can be put here.")
    ;

    po::options_description info_args("Information");
    bool help, version;
    info_args.add_options()
        ("help", po::bool_switch(&help)->default_value(false), "Display usage summary.")
        ("version", po::bool_switch(&version)->default_value(false), "Display program version.")
    ;

    po::options_description args;
    args.add(in_args).add(out_args).add(search_args).add(config_args).add(info_args);

    // Remain empty to prevent any positional option being provided 
    po::positional_options_description positional_args;

    try {
        // variable_map tells whether an option is explicity set in argv
        po::variables_map vm;
        po::store(
            po::command_line_parser(argc, argv)
                .options(args)
                .style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing)
                .positional(positional_args)
                .run(),
            vm);
        po::notify(vm);

        /*
         * Since required options are manually verified instead of using po::value.required(),
         * --version or --help can occur alone while required options are missing.
         */
        if (version) {
            std::cout << bd::VERSION_STR << std::endl;
            return 0;
        } else if (help) {
            std::cout << bd::VERSION_STR << '\n';
            std::cout << args << '\n';
            return 0;
        }

        // Read options set in config file
        if (vm.count("config") > 0) {
            po::options_description shared_args;
            shared_args.add(in_args).add(out_args).add(search_args);
            auto in = bd::open_for_read(config_file);
            po::store(po::parse_config_file(*in, shared_args), vm);
            po::notify(vm);
        }

        // Check inputs
        if (receptor_file.empty()) {
            throw po::required_option("receptor_file");
        }
        if (nthreads < 0) {
            throw po::invalid_option_value("nthreads");
        }

        // Check outputs
        if (output_dir.empty()) {
            throw po::invalid_option_value("output_dir");
        }
        if (log_file.empty()) {
            throw po::invalid_option_value("log_file");
        }
        if (verbosity < 0 || verbosity > 2) {
            throw po::invalid_option_value("verbosity");
        }

        // Check precalculation region
        if (vm.count("center_x") == 0) {
            throw po::required_option("center_x");
        }
        if (vm.count("center_y") == 0) {
            throw po::required_option("center_y");
        }
        if (vm.count("center_z") == 0) {
            throw po::required_option("center_z");
        }
        if (vm.count("size_x") == 0) {
            throw po::required_option("size_x");
        }
        if (vm.count("size_y") == 0) {
            throw po::required_option("size_y");
        }
        if (vm.count("size_z") == 0) {
            throw po::required_option("size_z");
        }
        if CHECK_PARAMETER_NON_POSITIVE(size_x) {
            throw po::invalid_option_value("size_x");
        }
        if CHECK_PARAMETER_NON_POSITIVE(size_y) {
            throw po::invalid_option_value("size_y");
        }
        if CHECK_PARAMETER_NON_POSITIVE(size_z) {
            throw po::invalid_option_value("size_z");
        }
        if CHECK_PARAMETER_NON_POSITIVE(grid_spacing - 0.01) {
            throw po::invalid_option_value("spacing");
        }
    } catch (po::error_with_option_name & e) {
        std::cerr << e.what() << std::endl;
        std::cout << args << std::endl;
        return 1;
    }

    // Let's generate receptor cache here
    try {
        bd::setup_global_logging(log_file, verbosity);
        bd::ReusableEngine engine(receptor_file, sf_file, nthreads);
        engine.set_box(center_x, center_y, center_z, size_x, size_y, size_z);
        engine.generate_maps(grid_spacing, max_size);
        engine.dump_maps(output_dir);
        return 0;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "An unknown error occurred." << std::endl;
        return 1;
    }
}
