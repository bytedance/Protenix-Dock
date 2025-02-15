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
    std::string sf_file, receptor_file, ligand_index_file;
    in_args.add_options()
        ("scoring_function", po::value<std::string>(&sf_file), "A Json file defining the scoring function. If empty, v6 is enabled.")
        ("receptor_file", po::value<std::string>(&receptor_file), "A Json file describing the receptor.")
        ("ligand_index", po::value<std::string>(&ligand_index_file), "Text file where each line is Json file describing ligand.")
    ;

    po::options_description out_args("Outputs");
    std::string output_dir, log_file;
    int verbosity;
    out_args.add_options()
        ("output_dir", po::value<std::string>(&output_dir)->default_value("docked"), "Output folder for optimized ligand poses.")
        ("log_file", po::value<std::string>(&log_file)->default_value("-"), "Text file for log messages. \"-\" is short for STDOUT.")
        ("verbosity", po::value<int>(&verbosity)->default_value(1), "Verbosity level. (0=warning, 1=info, 2=debug)")
    ;

    po::options_description search_args("Precalculation");
    double center_x, center_y, center_z;
    double size_x, size_y, size_z;
    double grid_spacing, slope;
    std::string map_dir;
    search_args.add_options()
        ("map_dir", po::value<std::string>(&map_dir), "Map folder for precalculated grids of receptor.")
        ("center_x", po::value<double>(&center_x), "X coordinate of box center (Angstrom).")
        ("center_y", po::value<double>(&center_y), "Y coordinate of box center (Angstrom).")
        ("center_z", po::value<double>(&center_z), "Z coordinate of box center (Angstrom).")
        ("size_x", po::value<double>(&size_x), "Box size in X dimension (Angstrom).")
        ("size_y", po::value<double>(&size_y), "Box size in Y dimension (Angstrom).")
        ("size_z", po::value<double>(&size_z), "Box size in Z dimension (Angstrom).")
        ("spacing", po::value<double>(&grid_spacing)->default_value(0.175), "Grid spacing. When value<0.01, disable it. (Angstrom)")
        ("penalty", po::value<double>(&slope)->default_value(1e6), "Penalty factor for ligand atoms out of box. (kJ/mol/Angstrom)")
    ;

    po::options_description algo_args("Algorithm");
    int seed;
    int exhaustiveness;
    int max_nsteps;
    int relax_nsteps;
    int num_modes;
    double min_rmsd;
    double mc_prune_energy_threshold;
    algo_args.add_options()
        ("seed", po::value<int>(&seed)->default_value(0), "Random seed for reproducibility.")
        ("exhaustiveness", po::value<int>(&exhaustiveness)->default_value(256), "Count of walkers for conformer space.")
        ("max_nsteps", po::value<int>(&max_nsteps)->default_value(200), "Max mutation count for one round of Monte-Carlo search.")
        ("relax_nsteps", po::value<int>(&relax_nsteps)->default_value(40), "Optimization steps to relax a mutated conformer.")
        ("num_modes", po::value<int>(&num_modes)->default_value(8), "Maximum binding poses to generate.")
        ("min_rmsd", po::value<double>(&min_rmsd)->default_value(0.5), "Minimum RMSD between generated poses.")
        ("mc_prune_energy_threshold", po::value<double>(&mc_prune_energy_threshold)->default_value(-1.), "MMenergy threshold for MC prune. When value<0, disable MC prune.")
    ;

    po::options_description device_args("Computing device");
    int num_readers, num_writers, nthreads;
    device_args.add_options()
        ("num_readers", po::value<int>(&num_readers)->default_value(1), "Thread count for parsing ligand files.")
        ("num_writers", po::value<int>(&num_writers)->default_value(1), "Thread count for writing results.")
        ("nthreads", po::value<int>(&nthreads)->default_value(0), "Thread count for docking tasks. 0 means all available cores.")
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
    args.add(in_args).add(out_args).add(search_args).add(algo_args).add(device_args).add(config_args).add(info_args);

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
            shared_args.add(in_args).add(out_args).add(search_args).add(algo_args);
            auto in = bd::open_for_read(config_file);
            po::store(po::parse_config_file(*in, shared_args), vm);
            po::notify(vm);
        }

        // Check inputs
        if (receptor_file.empty()) {
            throw po::required_option("receptor_file");
        }
        if (vm.count("ligand_index") == 0 || ligand_index_file.empty()) {
            throw po::required_option("ligand_index");
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
        if (map_dir.empty()) {
            /**
             * No matter whether receptor cache is enabled,
             * box information is always required for intial conformer genration.
             */
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
            if CHECK_PARAMETER_NON_POSITIVE(slope) {
                throw po::invalid_option_value("slope");
            }
        } else {
            // Box information is provided by map files
            if CHECK_PARAMETER_NON_POSITIVE(slope) {
                throw po::invalid_option_value("slope");
            }
        }

        // Check algorithm
        if (seed < 0) {
            throw po::invalid_option_value("seed");
        }
        if (exhaustiveness < 1) {
            throw po::invalid_option_value("exhaustiveness");
        }
        if (max_nsteps < 0) {
            throw po::invalid_option_value("max_nsteps");
        }
        if (relax_nsteps < 0) {
            throw po::invalid_option_value("relax_nsteps");
        }
        if (num_modes < 1) {
            throw po::invalid_option_value("num_modes");
        }
        if (min_rmsd < 0) {
            throw po::invalid_option_value("min_rmsd");
        }

        // Check computing device
        if (num_readers < 1) {
            throw po::invalid_option_value("num_readers");
        }
        if (num_writers < 1) {
            throw po::invalid_option_value("num_writers");
        }
        if (nthreads < 0) {
            throw po::invalid_option_value("nthreads");
        }
    } catch (po::error_with_option_name & e) {
        std::cerr << e.what() << std::endl;
        std::cout << args << std::endl;
        return 1;
    }

    // Let's dock multiple ligands here
    try {
        bd::setup_global_logging(log_file, verbosity);
        bd::ReusableEngine engine(receptor_file, sf_file, nthreads);
        if (!map_dir.empty()) {
            engine.load_maps(map_dir);
        } else {
            engine.set_box(center_x, center_y, center_z, size_x, size_y, size_z);
            if CHECK_PARAMETER_POSITIVE(grid_spacing) {
                engine.generate_maps(grid_spacing);
            }
        }
        engine.search(
            ligand_index_file,
            output_dir,
            seed,
            exhaustiveness,
            max_nsteps,
            relax_nsteps,
            num_modes,
            min_rmsd,
            mc_prune_energy_threshold,
            slope
        );
        return 0;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "An unknown error occurred." << std::endl;
        return 1;
    }
}
