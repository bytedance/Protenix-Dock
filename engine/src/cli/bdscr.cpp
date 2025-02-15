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
#include "bytedock/version.h"

namespace bd = bytedock;
namespace po = boost::program_options;

inline bool endswith(const std::string& value, const std::string& suffix) {
    if (value.size() < suffix.size()) return false;
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

int main(int argc, char* argv[]) {
    po::options_description in_args("Inputs");
    std::string sf_file, receptor_file, ligand_index_file;
    int nthreads;
    in_args.add_options()
        ("scoring_function", po::value<std::string>(&sf_file), "A Json file defining the scoring function. If empty, v6 is enabled.")
        ("receptor_file", po::value<std::string>(&receptor_file), "A Json file describing the receptor.")
        ("ligand_index", po::value<std::string>(&ligand_index_file), "Text file where each line is Json file describing ligand.")
        ("nthreads", po::value<int>(&nthreads)->default_value(0), "Thread count for evaluations. 0 means all available cores.")
    ;

    po::options_description out_args("Outputs");
    std::string output_file, log_file;
    int verbosity;
    out_args.add_options()
        ("output_file", po::value<std::string>(&output_file)->default_value("scores.csv"), "Output CSV file of ligand scores.")
        ("log_file", po::value<std::string>(&log_file)->default_value("-"), "Text file for log messages. \"-\" is short for STDOUT.")
        ("verbosity", po::value<int>(&verbosity)->default_value(1), "Verbosity level. (0=warning, 1=info, 2=debug)")
    ;

    po::options_description info_args("Information");
    bool help, version;
    info_args.add_options()
        ("help", po::bool_switch(&help)->default_value(false), "Display usage summary.")
        ("version", po::bool_switch(&version)->default_value(false), "Display program version.")
    ;

    po::options_description args;
    args.add(in_args).add(out_args).add(info_args);

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

        // Check inputs
        if (receptor_file.empty()) {
            throw po::required_option("receptor_file");
        }
        if (vm.count("ligand_index") == 0 || ligand_index_file.empty()) {
            throw po::required_option("ligand_index");
        }
        if (nthreads < 0) {
            throw po::invalid_option_value("nthreads");
        }

        // Check outputs
        if (!endswith(output_file, "csv")) {
            throw po::invalid_option_value("output_file");
        }
        if (log_file.empty()) {
            throw po::invalid_option_value("log_file");
        }
        if (verbosity < 0 || verbosity > 2) {
            throw po::invalid_option_value("verbosity");
        }
    } catch (po::error_with_option_name & e) {
        std::cerr << e.what() << std::endl;
        std::cout << args << std::endl;
        return 1;
    }

    // Let's evaluate scores of multiple ligands here
    try {
        bd::setup_global_logging(log_file, verbosity);
        bd::ReusableEngine engine(receptor_file, sf_file, nthreads);
        engine.evaluate(ligand_index_file, output_file);
        return 0;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "An unknown error occurred." << std::endl;
        return 1;
    }
}
