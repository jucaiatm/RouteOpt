#include "vurpm.hpp"

#include <bit>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

namespace RouteOpt::Application::VURPM {
namespace {

[[noreturn]] void throwGurobi(GRBmodel *model, int error, const std::string &context) {
    const char *message = nullptr;
    if (model != nullptr) {
        message = GRBgeterrormsg(GRBgetenv(model));
    }
    std::ostringstream oss;
    oss << context << " failed with Gurobi error " << error;
    if (message != nullptr) oss << ": " << message;
    throw std::runtime_error(oss.str());
}

void checkGurobi(GRBmodel *model, int error, const std::string &context) {
    if (error != 0) throwGurobi(model, error, context);
}

void checkGurobiEnv(GRBenv *environment, int error, const std::string &context) {
    if (error == 0) return;
    std::ostringstream oss;
    oss << context << " failed with Gurobi error " << error;
    if (environment != nullptr) {
        const char *message = GRBgeterrormsg(environment);
        if (message != nullptr) oss << ": " << message;
    }
    throw std::runtime_error(oss.str());
}

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string uppercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return value;
}

std::vector<std::string> tokenize(std::string line) {
    const auto hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    for (char &c : line) {
        if (c == ':' || c == ',' || c == ';' || c == '=') c = ' ';
    }
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) tokens.push_back(token);
    return tokens;
}

double euclidean(const Point &a, const Point &b) {
    return std::hypot(a.x - b.x, a.y - b.y);
}

std::uint64_t fullMask(int n) {
    if (n <= 0 || n > 63) throw std::runtime_error("VURP-M supports 1..63 platforms");
    return n == 64 ? std::numeric_limits<std::uint64_t>::max() : ((std::uint64_t{1} << n) - 1U);
}

int popcount(std::uint64_t mask) {
    return std::popcount(mask);
}

int leastBitIndex(std::uint64_t mask) {
    return static_cast<int>(std::countr_zero(mask));
}

std::string gurobiName(const std::string &name) {
    if (name.size() <= 180) return name;
    const std::size_t hash = std::hash<std::string>{}(name);
    return name.substr(0, 140) + "_" + std::to_string(hash);
}

std::string quantizedDirection(double x, double y) {
    constexpr double scale = 1e9;
    const auto qx = static_cast<long long>(std::llround(x * scale));
    const auto qy = static_cast<long long>(std::llround(y * scale));
    return std::to_string(qx) + ":" + std::to_string(qy);
}

std::uint64_t slotPatternKey(int slot, int pattern_id) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(slot)) << 32U) |
           static_cast<std::uint32_t>(pattern_id);
}

bool fractional(double value) {
    return std::abs(value - std::round(value)) > kIntegralityTolerance;
}

struct CandidateState {
    double reduced_cost{};
    std::uint64_t mask{};
    int end{-1};
    int start{-1};
    double distance{};
};

struct WorseCandidate {
    bool operator()(const CandidateState &a, const CandidateState &b) const {
        return a.reduced_cost < b.reduced_cost;
    }
};

int addRawVar(GRBmodel *model, double obj, double lb, double ub, char type, const std::string &name) {
    int num_vars = 0;
    checkGurobi(model, GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &num_vars), "GRBgetintattr(NumVars)");
    checkGurobi(model, GRBaddvar(model, 0, nullptr, nullptr, obj, lb, ub, type, name.c_str()),
                 "GRBaddvar(" + name + ")");
    checkGurobi(model, GRBupdatemodel(model), "GRBupdatemodel after variable");
    return num_vars;
}

void addRawLinearConstraint(GRBmodel *model,
                            const std::vector<std::pair<int, double>> &terms,
                            char sense,
                            double rhs,
                            const std::string &name) {
    std::vector<int> indices;
    std::vector<double> values;
    indices.reserve(terms.size());
    values.reserve(terms.size());
    for (const auto &[index, coefficient] : terms) {
        if (std::abs(coefficient) <= 1e-14) continue;
        indices.push_back(index);
        values.push_back(coefficient);
    }
    checkGurobi(model,
                 GRBaddconstr(model, static_cast<int>(indices.size()),
                              indices.empty() ? nullptr : indices.data(),
                              values.empty() ? nullptr : values.data(),
                              sense, rhs, name.c_str()),
                 "GRBaddconstr(" + name + ")");
}

void addNormVarConst(GRBmodel *model,
                     int x_var, int y_var, const Point &constant,
                     int time_var, double speed,
                     const std::string &name) {
    std::array<int, 2> lind{x_var, y_var};
    std::array<double, 2> lval{-2.0 * constant.x, -2.0 * constant.y};
    std::array<int, 3> qrow{x_var, y_var, time_var};
    std::array<int, 3> qcol{x_var, y_var, time_var};
    std::array<double, 3> qval{1.0, 1.0, -speed * speed};
    const double rhs = -(constant.x * constant.x + constant.y * constant.y);
    checkGurobi(model,
                 GRBaddqconstr(model,
                               static_cast<int>(lind.size()), lind.data(), lval.data(),
                               static_cast<int>(qrow.size()), qrow.data(), qcol.data(), qval.data(),
                               GRB_LESS_EQUAL, rhs, name.c_str()),
                 "GRBaddqconstr(" + name + ")");
}

void addNormVarVar(GRBmodel *model,
                   int x1, int y1, int x2, int y2,
                   int time_var, double speed,
                   const std::string &name) {
    std::array<int, 7> qrow{x1, y1, x2, y2, x1, y1, time_var};
    std::array<int, 7> qcol{x1, y1, x2, y2, x2, y2, time_var};
    std::array<double, 7> qval{1.0, 1.0, 1.0, 1.0, -2.0, -2.0, -speed * speed};
    checkGurobi(model,
                 GRBaddqconstr(model, 0, nullptr, nullptr,
                               static_cast<int>(qrow.size()), qrow.data(), qcol.data(), qval.data(),
                               GRB_LESS_EQUAL, 0.0, name.c_str()),
                 "GRBaddqconstr(" + name + ")");
}

} // namespace

bool OmegaBox::contains(const Point &p, double tol) const {
    return p.x >= xmin - tol && p.x <= xmax + tol && p.y >= ymin - tol && p.y <= ymax + tol;
}

Instance Instance::read(const std::string &path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open VURP-M instance: " + path);

    Instance instance;
    bool reading_platforms = false;
    std::unordered_map<int, Point> by_id;
    std::string line;
    int line_number = 0;

    while (std::getline(input, line)) {
        ++line_number;
        auto tokens = tokenize(line);
        if (tokens.empty()) continue;
        const std::string key = uppercase(tokens.front());

        if (key == "PLATFORMS" || key == "NODE_COORD_SECTION") {
            reading_platforms = true;
            continue;
        }
        if (key == "EOF") break;

        if (reading_platforms) {
            if (tokens.size() < 3) {
                throw std::runtime_error("Invalid platform row at line " + std::to_string(line_number));
            }
            const int id = std::stoi(tokens[0]);
            by_id[id] = Point{std::stod(tokens[1]), std::stod(tokens[2])};
            continue;
        }

        if (key == "NAME" && tokens.size() >= 2) {
            instance.name = tokens[1];
        } else if ((key == "DIMENSION" || key == "N") && tokens.size() >= 2) {
            instance.n = std::stoi(tokens[1]);
        } else if ((key == "VESSEL_SPEED" || key == "VV") && tokens.size() >= 2) {
            instance.vessel_speed = std::stod(tokens[1]);
        } else if ((key == "UAV_SPEED" || key == "VD") && tokens.size() >= 2) {
            instance.uav_speed = std::stod(tokens[1]);
        } else if ((key == "ENDURANCE" || key == "B") && tokens.size() >= 2) {
            instance.endurance = std::stod(tokens[1]);
        } else if (key == "DEPOT" && tokens.size() >= 3) {
            instance.depot = Point{std::stod(tokens[1]), std::stod(tokens[2])};
        } else if (key == "OMEGA" && tokens.size() >= 5) {
            instance.omega = OmegaBox{std::stod(tokens[1]), std::stod(tokens[2]),
                                      std::stod(tokens[3]), std::stod(tokens[4])};
            instance.omega_explicit = true;
        } else {
            throw std::runtime_error("Unknown instance key at line " + std::to_string(line_number) + ": " + tokens[0]);
        }
    }

    if (instance.n == 0) instance.n = static_cast<int>(by_id.size());
    instance.platform.resize(instance.n);
    for (int i = 1; i <= instance.n; ++i) {
        const auto it = by_id.find(i);
        if (it == by_id.end()) {
            throw std::runtime_error("Missing platform id " + std::to_string(i));
        }
        instance.platform[i - 1] = it->second;
    }

    instance.preprocess();
    instance.validate();
    return instance;
}

void Instance::preprocess() {
    if (!omega_explicit && !platform.empty()) {
        double xmin = depot.x;
        double xmax = depot.x;
        double ymin = depot.y;
        double ymax = depot.y;
        for (const auto &p : platform) {
            xmin = std::min(xmin, p.x);
            xmax = std::max(xmax, p.x);
            ymin = std::min(ymin, p.y);
            ymax = std::max(ymax, p.y);
        }
        const double margin = std::max(1.0, uav_speed * endurance);
        omega = OmegaBox{xmin - margin, xmax + margin, ymin - margin, ymax + margin};
    }

    distance.assign(n, std::vector<double>(n, 0.0));
    depot_distance.assign(n, 0.0);
    conflict_mask.assign(n, 0U);
    for (int i = 0; i < n; ++i) {
        depot_distance[i] = euclidean(depot, platform[i]);
        for (int j = 0; j < n; ++j) {
            distance[i][j] = euclidean(platform[i], platform[j]);
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (2.0 * distance[i][j] > (uav_speed + vessel_speed) * endurance + kTolerance) {
                conflict_mask[i] |= std::uint64_t{1} << j;
                conflict_mask[j] |= std::uint64_t{1} << i;
            }
        }
    }

    const std::array<Point, 4> corners{{
        {omega.xmin, omega.ymin}, {omega.xmin, omega.ymax},
        {omega.xmax, omega.ymin}, {omega.xmax, omega.ymax}
    }};
    big_m_distance = 0.0;
    for (const auto &corner : corners) {
        big_m_distance = std::max(big_m_distance, std::hypot(corner.x, corner.y));
        for (const auto &p : platform) {
            big_m_distance = std::max(big_m_distance, euclidean(corner, p));
        }
    }
    big_m_distance = std::max(1.0, big_m_distance + 1.0);
}

void Instance::validate() const {
    if (n <= 0 || n > 63) throw std::runtime_error("DIMENSION must be in 1..63");
    if (static_cast<int>(platform.size()) != n) throw std::runtime_error("Platform count does not match DIMENSION");
    if (vessel_speed <= 0.0) throw std::runtime_error("VESSEL_SPEED must be positive");
    if (uav_speed <= 0.0) throw std::runtime_error("UAV_SPEED must be positive");
    if (endurance <= 0.0) throw std::runtime_error("ENDURANCE must be positive");
    if (omega.xmin > omega.xmax || omega.ymin > omega.ymax) throw std::runtime_error("Invalid OMEGA box");
    if (!omega.contains(depot)) throw std::runtime_error("DEPOT must lie in OMEGA");
    for (int i = 0; i < n; ++i) {
        if (!omega.contains(platform[i])) {
            throw std::runtime_error("Platform " + std::to_string(i + 1) + " must lie in OMEGA");
        }
    }
}

std::string RoutePattern::key() const {
    std::ostringstream oss;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        if (i != 0) oss << '-';
        oss << sequence[i];
    }
    return oss.str();
}

int PatternPool::add(const RoutePattern &pattern) {
    const std::string pattern_key = pattern.key();
    const auto existing = index_.find(pattern_key);
    if (existing != index_.end()) return existing->second;
    RoutePattern copy = pattern;
    copy.id = static_cast<int>(patterns_.size());
    patterns_.push_back(copy);
    index_.emplace(pattern_key, copy.id);
    return copy.id;
}

const RoutePattern &PatternPool::get(int id) const {
    if (id < 0 || id >= static_cast<int>(patterns_.size())) throw std::out_of_range("Invalid route pattern id");
    return patterns_[id];
}

int PatternPool::size() const {
    return static_cast<int>(patterns_.size());
}

std::optional<int> PatternPool::find(const RoutePattern &pattern) const {
    const auto it = index_.find(pattern.key());
    if (it == index_.end()) return std::nullopt;
    return it->second;
}

ExactPricing::ExactPricing(const Instance &instance, int max_exact_customers)
    : instance_(instance), max_exact_customers_(max_exact_customers) {
}

bool ExactPricing::subsetConflictFree(std::uint64_t mask) const {
