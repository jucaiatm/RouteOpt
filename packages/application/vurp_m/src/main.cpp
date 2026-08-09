#include "vurpm.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using RouteOpt::Application::VURPM::BpcParameters;
using RouteOpt::Application::VURPM::BpcSolver;
using RouteOpt::Application::VURPM::Instance;

namespace {

void printUsage(const char *program) {
    std::cerr << "Usage: " << program << " <instance> [options]\n"
              << "Options:\n"
              << "  --time-limit <seconds>\n"
              << "  --threads <count>\n"
              << "  --pricing-columns <count>\n"
              << "  --oa-cuts <count>\n"
              << "  --max-exact-customers <count>\n"
              << "  --verbose-lp\n"
              << "  --verbose-socp\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    try {
        BpcParameters parameters;
        const std::string instance_path = argv[1];
        for (int i = 2; i < argc; ++i) {
            const std::string option = argv[i];
            auto requireValue = [&](const std::string &name) -> std::string {
                if (i + 1 >= argc) throw std::runtime_error("Missing value for " + name);
                return argv[++i];
            };
            if (option == "--time-limit") {
                parameters.time_limit = std::stod(requireValue(option));
            } else if (option == "--threads") {
                parameters.threads = std::stoi(requireValue(option));
            } else if (option == "--pricing-columns") {
                parameters.max_pricing_columns_per_slot = std::stoi(requireValue(option));
            } else if (option == "--oa-cuts") {
                parameters.max_oa_cuts_per_round = std::stoi(requireValue(option));
            } else if (option == "--max-exact-customers") {
                parameters.max_exact_pricing_customers = std::stoi(requireValue(option));
            } else if (option == "--verbose-lp") {
                parameters.verbose_lp = true;
            } else if (option == "--verbose-socp") {
                parameters.verbose_socp = true;
            } else if (option == "--help" || option == "-h") {
                printUsage(argv[0]);
                return EXIT_SUCCESS;
            } else {
                throw std::runtime_error("Unknown option: " + option);
            }
        }

        Instance instance = Instance::read(instance_path);
        BpcSolver solver(std::move(instance), parameters);
        solver.solve();
        solver.printResult();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "VURP-M BPC error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
