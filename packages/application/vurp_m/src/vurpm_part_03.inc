                if (fixing.value == 1) {
                    if (result.forced_first >= 0 && result.forced_first != fixing.platform) result.impossible = true;
                    result.forced_first = fixing.platform;
                    result.required_mask |= bit;
                } else {
                    result.forbidden_first_mask |= bit;
                }
                break;
            case BranchVarKind::LastPlatform:
                if (fixing.value == 1) {
                    if (result.forced_last >= 0 && result.forced_last != fixing.platform) result.impossible = true;
                    result.forced_last = fixing.platform;
                    result.required_mask |= bit;
                } else {
                    result.forbidden_last_mask |= bit;
                }
                break;
        }
    }
    if ((result.required_mask & result.forbidden_mask) != 0U) result.impossible = true;
    if (result.forced_first >= 0 && ((result.forbidden_first_mask >> result.forced_first) & 1U)) result.impossible = true;
    if (result.forced_last >= 0 && ((result.forbidden_last_mask >> result.forced_last) & 1U)) result.impossible = true;
    return result;
}

std::vector<LinearCut> MasterProblem::separateOuterApproximation(
    const MasterSolution &solution,
    int maximum_cuts) const {

    std::vector<LinearCut> cuts;
    if (!solution.optimal || solution.phase_one || maximum_cuts <= 0) return cuts;
    const auto &x = solution.values;
    auto consider = [&](int type, int slot, double dx, double dy, double rhs) {
        if (static_cast<int>(cuts.size()) >= maximum_cuts) return;
        const double norm = std::hypot(dx, dy);
        if (norm <= rhs + kOATolerance || norm <= 1e-14) return;
        cuts.push_back(makeOACut(type, slot, Point{dx / norm, dy / norm}));
    };

    consider(0, 0,
             x[layout_.wx[0]] - instance_.depot.x,
             x[layout_.wy[0]] - instance_.depot.y,
             instance_.vessel_speed * x[layout_.ts]);

    const int last_slot = layout_.slots - 1;
    consider(1, last_slot,
             x[layout_.qx[last_slot]] - instance_.depot.x,
             x[layout_.qy[last_slot]] - instance_.depot.y,
             instance_.vessel_speed * x[layout_.tf]);

    for (int k = 0; k < layout_.slots; ++k) {
        consider(2, k,
                 x[layout_.qx[k]] - x[layout_.wx[k]],
                 x[layout_.qy[k]] - x[layout_.wy[k]],
                 instance_.vessel_speed * x[layout_.td[k]]);
        if (k + 1 < layout_.slots) {
            consider(3, k,
                     x[layout_.wx[k + 1]] - x[layout_.qx[k]],
                     x[layout_.wy[k + 1]] - x[layout_.qy[k]],
                     instance_.vessel_speed * x[layout_.tb[k]]);
        }

        Point first_point{};
        Point last_point{};
        for (int i = 0; i < layout_.n; ++i) {
            const double first_value = x[layout_.at(layout_.first, k, i)];
            const double last_value = x[layout_.at(layout_.last, k, i)];
            first_point.x += instance_.platform[i].x * first_value;
            first_point.y += instance_.platform[i].y * first_value;
            last_point.x += instance_.platform[i].x * last_value;
            last_point.y += instance_.platform[i].y * last_value;
        }
        const double relaxed_first_rhs = x[layout_.alpha[k]] + instance_.big_m_distance * (1.0 - x[layout_.z[k]]);
        const double relaxed_last_rhs = x[layout_.beta[k]] + instance_.big_m_distance * (1.0 - x[layout_.z[k]]);
        consider(4, k,
                 x[layout_.wx[k]] - first_point.x,
                 x[layout_.wy[k]] - first_point.y,
                 relaxed_first_rhs);
        consider(5, k,
                 last_point.x - x[layout_.qx[k]],
                 last_point.y - x[layout_.qy[k]],
                 relaxed_last_rhs);
    }
    return cuts;
}

bool MasterProblem::staticVariablesIntegral(const MasterSolution &solution) const {
    if (!solution.optimal || solution.phase_one) return false;
    for (const int index : layout_.z) if (fractional(solution.values[index])) return false;
    for (const int index : layout_.y) if (fractional(solution.values[index])) return false;
    for (const int index : layout_.first) if (fractional(solution.values[index])) return false;
    for (const int index : layout_.last) if (fractional(solution.values[index])) return false;
    return true;
}

std::optional<BranchChoice> MasterProblem::chooseBranch(const MasterSolution &solution) const {
    if (!solution.optimal || solution.phase_one) return std::nullopt;
    std::optional<BranchChoice> best;
    auto inspect = [&](BranchVarKind kind, int slot, int platform, int column) {
        const double value = solution.values[column];
        if (!fractional(value)) return;
        const double score = std::abs(value - 0.5);
        if (!best || score < best->fractionality) {
            best = BranchChoice{Fixing{kind, slot, platform, 0}, Fixing{kind, slot, platform, 1}, score};
        }
    };

    for (int k = 0; k < layout_.slots; ++k) inspect(BranchVarKind::SortieActive, k, -1, layout_.z[k]);
    for (int k = 0; k < layout_.slots; ++k) {
        for (int i = 0; i < layout_.n; ++i) {
            inspect(BranchVarKind::PlatformAssignment, k, i, layout_.at(layout_.y, k, i));
        }
    }
    for (int k = 0; k < layout_.slots; ++k) {
        for (int i = 0; i < layout_.n; ++i) {
            inspect(BranchVarKind::FirstPlatform, k, i, layout_.at(layout_.first, k, i));
            inspect(BranchVarKind::LastPlatform, k, i, layout_.at(layout_.last, k, i));
        }
    }
    return best;
}

std::string MasterProblem::partitionKey(const MasterSolution &solution) const {
    std::string key;
    key.reserve(layout_.z.size() + layout_.y.size() + layout_.first.size() + layout_.last.size() + 8);
    auto append = [&](const std::vector<int> &indices) {
        for (const int index : indices) key.push_back(solution.values[index] >= 0.5 ? '1' : '0');
        key.push_back('|');
    };
    append(layout_.z);
    append(layout_.y);
    append(layout_.first);
    append(layout_.last);
    return key;
}

SocpOracle::SocpOracle(GRBenv *environment, const Instance &instance, bool verbose)
    : environment_(environment), instance_(instance), verbose_(verbose) {
}

SocpSolution SocpOracle::solve(const PartitionPlan &partition) const {
    SocpSolution solution;
    const int m = static_cast<int>(partition.sorties.size());
    if (m <= 0) {
        solution.proven_infeasible = true;
        return solution;
    }

    GRBmodel *model = nullptr;
    checkGurobiEnv(environment_, GRBnewmodel(environment_, &model, "VURP_M_SOCP", 0,
                                             nullptr, nullptr, nullptr, nullptr, nullptr),
                    "GRBnewmodel(VURP_M_SOCP)");
    try {
        checkGurobi(model, GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_OUTPUTFLAG, verbose_ ? 1 : 0),
                     "set SOCP OutputFlag");
        checkGurobi(model, GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_NONCONVEX, 0),
                     "set SOCP NonConvex");
        checkGurobi(model, GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE), "set SOCP sense");

        std::vector<int> wx(m), wy(m), qx(m), qy(m), alpha(m), beta(m), td(m);
        std::vector<int> tb(std::max(0, m - 1));
        for (int k = 0; k < m; ++k) {
            wx[k] = addRawVar(model, 0.0, instance_.omega.xmin, instance_.omega.xmax, GRB_CONTINUOUS,
                              "wx_" + std::to_string(k));
            wy[k] = addRawVar(model, 0.0, instance_.omega.ymin, instance_.omega.ymax, GRB_CONTINUOUS,
                              "wy_" + std::to_string(k));
            qx[k] = addRawVar(model, 0.0, instance_.omega.xmin, instance_.omega.xmax, GRB_CONTINUOUS,
                              "qx_" + std::to_string(k));
            qy[k] = addRawVar(model, 0.0, instance_.omega.ymin, instance_.omega.ymax, GRB_CONTINUOUS,
                              "qy_" + std::to_string(k));
            alpha[k] = addRawVar(model, 0.0, 0.0, instance_.max_internal_distance(), GRB_CONTINUOUS,
                                 "alpha_" + std::to_string(k));
            beta[k] = addRawVar(model, 0.0, 0.0, instance_.max_internal_distance(), GRB_CONTINUOUS,
                                "beta_" + std::to_string(k));
            td[k] = addRawVar(model, 1.0, 0.0, instance_.endurance, GRB_CONTINUOUS,
                              "td_" + std::to_string(k));
            if (k + 1 < m) {
                tb[k] = addRawVar(model, 1.0, 0.0, kInfinity, GRB_CONTINUOUS,
                                  "tb_" + std::to_string(k));
            }
        }
        const int ts = addRawVar(model, 1.0, 0.0, kInfinity, GRB_CONTINUOUS, "Ts");
        const int tf = addRawVar(model, 1.0, 0.0, kInfinity, GRB_CONTINUOUS, "Tf");

        for (int k = 0; k < m; ++k) {
            addRawLinearConstraint(model,
                                   {{alpha[k], 1.0}, {beta[k], 1.0}, {td[k], -instance_.uav_speed}},
                                   GRB_LESS_EQUAL, -partition.sorties[k].route.internal_distance,
                                   "uav_duration_" + std::to_string(k));
            addNormVarVar(model, qx[k], qy[k], wx[k], wy[k], td[k], instance_.vessel_speed,
                          "sync_" + std::to_string(k));
            const Point first_point = instance_.platform[partition.sorties[k].route.first];
            const Point last_point = instance_.platform[partition.sorties[k].route.last];
            addNormVarConst(model, wx[k], wy[k], first_point, alpha[k], 1.0,
                            "first_leg_" + std::to_string(k));
            addNormVarConst(model, qx[k], qy[k], last_point, beta[k], 1.0,
                            "last_leg_" + std::to_string(k));
            if (k + 1 < m) {
                addNormVarVar(model, wx[k + 1], wy[k + 1], qx[k], qy[k], tb[k], instance_.vessel_speed,
                              "bridge_" + std::to_string(k));
            }
        }
        addNormVarConst(model, wx[0], wy[0], instance_.depot, ts, instance_.vessel_speed, "depot_start");
        addNormVarConst(model, qx[m - 1], qy[m - 1], instance_.depot, tf, instance_.vessel_speed, "depot_finish");

        checkGurobi(model, GRBupdatemodel(model), "update SOCP");
        checkGurobi(model, GRBoptimize(model), "optimize SOCP");
        int status = 0;
        checkGurobi(model, GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status), "get SOCP status");
        if (status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD) {
            solution.proven_infeasible = true;
            GRBfreemodel(model);
            return solution;
        }
        if (status != GRB_OPTIMAL) {
            throw std::runtime_error("SOCP oracle did not reach optimality; status=" + std::to_string(status));
        }

        solution.feasible = true;
        checkGurobi(model, GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &solution.objective), "get SOCP objective");
        int num_vars = 0;
        checkGurobi(model, GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &num_vars), "get SOCP NumVars");
        std::vector<double> x(num_vars);
        checkGurobi(model, GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, num_vars, x.data()), "get SOCP X");
        solution.launch.resize(m);
        solution.recovery.resize(m);
        solution.sortie_time.resize(m);
        solution.bridge_time.resize(std::max(0, m - 1));
        for (int k = 0; k < m; ++k) {
            solution.launch[k] = Point{x[wx[k]], x[wy[k]]};
            solution.recovery[k] = Point{x[qx[k]], x[qy[k]]};
            solution.sortie_time[k] = x[td[k]];
            if (k + 1 < m) solution.bridge_time[k] = x[tb[k]];
        }
        solution.start_time = x[ts];
        solution.finish_time = x[tf];
        GRBfreemodel(model);
        return solution;
    } catch (...) {
        if (model != nullptr) GRBfreemodel(model);
        throw;
    }
}

BpcSolver::BpcSolver(Instance instance, BpcParameters parameters)
    : instance_(std::move(instance)),
      parameters_(parameters),
      pricing_(instance_, parameters_.max_exact_pricing_customers) {
    instance_.validate();
    initializeEnvironment();
    socp_oracle_ = std::make_unique<SocpOracle>(environment_, instance_, parameters_.verbose_socp);
    slot_patterns_.resize(instance_.n);
}

BpcSolver::~BpcSolver() {
    socp_oracle_.reset();
    if (environment_ != nullptr) GRBfreeenv(environment_);
}

void BpcSolver::initializeEnvironment() {
    int error = GRBloadenv(&environment_, nullptr);
    checkGurobiEnv(environment_, error, "GRBloadenv");
    checkGurobiEnv(environment_, GRBsetintparam(environment_, GRB_INT_PAR_OUTPUTFLAG, 0), "set global OutputFlag");
    checkGurobiEnv(environment_, GRBsetintparam(environment_, GRB_INT_PAR_THREADS, parameters_.threads),
                    "set global Threads");
}

void BpcSolver::initializePatterns() {
    for (int i = 0; i < instance_.n; ++i) {
        RoutePattern singleton;
        singleton.sequence = {i};
        singleton.visited_mask = std::uint64_t{1} << i;
        singleton.first = i;
        singleton.last = i;
        singleton.internal_distance = 0.0;
        const int id = patterns_.add(singleton);
        for (int k = 0; k < instance_.n; ++k) slot_patterns_[k].push_back(id);
    }
    statistics_.generated_patterns = patterns_.size();
    statistics_.generated_slot_columns = instance_.n * instance_.n;
}

void BpcSolver::initializeIncumbent() {
    PartitionPlan partition;
    partition.z.assign(instance_.n, 1);
    partition.y.assign(instance_.n * instance_.n, 0);
    partition.first.assign(instance_.n * instance_.n, 0);
    partition.last.assign(instance_.n * instance_.n, 0);

    std::vector<int> unvisited(instance_.n);
    std::iota(unvisited.begin(), unvisited.end(), 0);
    Point current = instance_.depot;
    for (int k = 0; k < instance_.n; ++k) {
        auto best_it = unvisited.begin();
        double best_distance = kInfinity;
        for (auto it = unvisited.begin(); it != unvisited.end(); ++it) {
            const double distance = euclidean(current, instance_.platform[*it]);
            if (distance < best_distance) {
                best_distance = distance;
                best_it = it;
            }
        }
        const int customer = *best_it;
        unvisited.erase(best_it);
        current = instance_.platform[customer];

        RoutePattern route;
        route.sequence = {customer};
        route.visited_mask = std::uint64_t{1} << customer;
        route.first = customer;
        route.last = customer;
        route.internal_distance = 0.0;
        partition.sorties.push_back(SortiePlan{k, route.visited_mask, customer, customer, route});
        partition.y[k * instance_.n + customer] = 1;
        partition.first[k * instance_.n + customer] = 1;
        partition.last[k * instance_.n + customer] = 1;
    }

    ++statistics_.socp_calls;
    const SocpSolution socp = socp_oracle_->solve(partition);
    if (socp.feasible) {
        incumbent_.available = true;
        incumbent_.objective = socp.objective;
        incumbent_.partition = partition;
        incumbent_.continuous = socp;
        std::cout << "Initial certified UB: " << std::setprecision(12) << incumbent_.objective << '\n';
    }
}

bool BpcSolver::timeLimitReached() const {
    return elapsed() >= parameters_.time_limit;
}

double BpcSolver::elapsed() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time_).count();
}

bool BpcSolver::fixingSetContradictory(const std::vector<Fixing> &fixings) const {
    std::unordered_map<std::string, int> values;
    std::vector<int> forced_first(instance_.n, -1);
    std::vector<int> forced_last(instance_.n, -1);
    std::vector<int> fixed_z(instance_.n, -1);
    std::vector<std::uint64_t> required(instance_.n, 0U);
    std::vector<std::uint64_t> forbidden(instance_.n, 0U);

    for (const auto &fixing : fixings) {
        const std::string key = std::to_string(static_cast<int>(fixing.kind)) + ":" +
                                std::to_string(fixing.slot) + ":" + std::to_string(fixing.platform);
        const auto it = values.find(key);
        if (it != values.end() && it->second != fixing.value) return true;
        values[key] = fixing.value;
        if (fixing.kind == BranchVarKind::SortieActive) fixed_z[fixing.slot] = fixing.value;
        if (fixing.kind == BranchVarKind::PlatformAssignment) {
            if (fixing.value == 1) required[fixing.slot] |= std::uint64_t{1} << fixing.platform;
            else forbidden[fixing.slot] |= std::uint64_t{1} << fixing.platform;
        }
        if (fixing.kind == BranchVarKind::FirstPlatform && fixing.value == 1) {
            if (forced_first[fixing.slot] >= 0 && forced_first[fixing.slot] != fixing.platform) return true;
            forced_first[fixing.slot] = fixing.platform;
            required[fixing.slot] |= std::uint64_t{1} << fixing.platform;
