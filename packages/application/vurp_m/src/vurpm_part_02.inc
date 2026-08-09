                                    "wy_" + std::to_string(k));
        layout_.qx[k] = addVariable(0.0, instance_.omega.xmin, instance_.omega.xmax, GRB_CONTINUOUS,
                                    "qx_" + std::to_string(k));
        layout_.qy[k] = addVariable(0.0, instance_.omega.ymin, instance_.omega.ymax, GRB_CONTINUOUS,
                                    "qy_" + std::to_string(k));
        layout_.alpha[k] = addVariable(0.0, 0.0, instance_.max_internal_distance(), GRB_CONTINUOUS,
                                       "alpha_" + std::to_string(k));
        layout_.beta[k] = addVariable(0.0, 0.0, instance_.max_internal_distance(), GRB_CONTINUOUS,
                                      "beta_" + std::to_string(k));
        layout_.td[k] = addVariable(0.0, 0.0, instance_.endurance, GRB_CONTINUOUS,
                                    "td_" + std::to_string(k));
        if (k + 1 < k_count) {
            layout_.tb[k] = addVariable(0.0, 0.0, kInfinity, GRB_CONTINUOUS,
                                        "tb_" + std::to_string(k));
        }
    }
    layout_.ts = addVariable(0.0, 0.0, kInfinity, GRB_CONTINUOUS, "Ts");
    layout_.tf = addVariable(0.0, 0.0, kInfinity, GRB_CONTINUOUS, "Tf");
    layout_.cmax = addVariable(0.0, 0.0, kInfinity, GRB_CONTINUOUS, "Cmax");
    layout_.static_variable_count = next_column_index_;
}

void MasterProblem::buildBaseRows() {
    const int n = layout_.n;
    const int k_count = layout_.slots;

    std::vector<std::pair<int, double>> objective_balance{{layout_.cmax, 1.0}, {layout_.ts, -1.0},
                                                          {layout_.tf, -1.0}};
    for (const int td : layout_.td) objective_balance.emplace_back(td, -1.0);
    for (const int tb : layout_.tb) objective_balance.emplace_back(tb, -1.0);
    addConstraint(objective_balance, GRB_EQUAL, 0.0, "total_completion_time");

    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<int, double>> terms;
        for (int k = 0; k < k_count; ++k) terms.emplace_back(layout_.at(layout_.y, k, i), 1.0);
        addConstraint(terms, GRB_EQUAL, 1.0, "cover_" + std::to_string(i));
    }

    for (int k = 0; k < k_count; ++k) {
        std::vector<std::pair<int, double>> first_sum{{layout_.z[k], -1.0}};
        std::vector<std::pair<int, double>> last_sum{{layout_.z[k], -1.0}};
        for (int i = 0; i < n; ++i) {
            const int y = layout_.at(layout_.y, k, i);
            const int first = layout_.at(layout_.first, k, i);
            const int last = layout_.at(layout_.last, k, i);
            addConstraint({{y, 1.0}, {layout_.z[k], -1.0}}, GRB_LESS_EQUAL, 0.0,
                          "activate_" + std::to_string(i) + "_" + std::to_string(k));
            addConstraint({{first, 1.0}, {y, -1.0}}, GRB_LESS_EQUAL, 0.0,
                          "first_visit_" + std::to_string(i) + "_" + std::to_string(k));
            addConstraint({{last, 1.0}, {y, -1.0}}, GRB_LESS_EQUAL, 0.0,
                          "last_visit_" + std::to_string(i) + "_" + std::to_string(k));
            first_sum.emplace_back(first, 1.0);
            last_sum.emplace_back(last, 1.0);
        }
        addConstraint(first_sum, GRB_EQUAL, 0.0, "one_first_" + std::to_string(k));
        addConstraint(last_sum, GRB_EQUAL, 0.0, "one_last_" + std::to_string(k));
    }

    for (int k = 0; k + 1 < k_count; ++k) {
        addConstraint({{layout_.z[k + 1], 1.0}, {layout_.z[k], -1.0}}, GRB_LESS_EQUAL, 0.0,
                      "slot_order_" + std::to_string(k));
    }

    rows_.convexity.resize(k_count);
    rows_.visit_link.resize(k_count * n);
    rows_.first_link.resize(k_count * n);
    rows_.last_link.resize(k_count * n);
    rows_.length_link.resize(k_count);

    for (int k = 0; k < k_count; ++k) {
        rows_.convexity[k] = addConstraint({{layout_.z[k], -1.0}}, GRB_EQUAL, 0.0,
                                           "convexity_" + std::to_string(k));
        phase_one_rows_.emplace_back(rows_.convexity[k], "convexity_" + std::to_string(k));
        for (int i = 0; i < n; ++i) {
            const int pos = k * n + i;
            rows_.visit_link[pos] = addConstraint({{layout_.y[pos], -1.0}}, GRB_EQUAL, 0.0,
                                                   "visit_link_" + std::to_string(i) + "_" + std::to_string(k));
            rows_.first_link[pos] = addConstraint({{layout_.first[pos], -1.0}}, GRB_EQUAL, 0.0,
                                                   "first_link_" + std::to_string(i) + "_" + std::to_string(k));
            rows_.last_link[pos] = addConstraint({{layout_.last[pos], -1.0}}, GRB_EQUAL, 0.0,
                                                  "last_link_" + std::to_string(i) + "_" + std::to_string(k));
            phase_one_rows_.emplace_back(rows_.visit_link[pos], "visit_link_" + std::to_string(pos));
            phase_one_rows_.emplace_back(rows_.first_link[pos], "first_link_" + std::to_string(pos));
            phase_one_rows_.emplace_back(rows_.last_link[pos], "last_link_" + std::to_string(pos));
        }
        // The inequality is equivalent at an optimum because ell only makes the
        // duration constraints harder.  Its nonnegative dual is essential for
        // exact shortest-Hamilton-path pricing.
        rows_.length_link[k] = addConstraint({{layout_.ell[k], 1.0}}, GRB_GREATER_EQUAL, 0.0,
                                             "length_link_" + std::to_string(k));

        addConstraint({{layout_.alpha[k], 1.0}, {layout_.ell[k], 1.0}, {layout_.beta[k], 1.0},
                       {layout_.td[k], -instance_.uav_speed}},
                      GRB_LESS_EQUAL, 0.0, "uav_duration_" + std::to_string(k));
        addConstraint({{layout_.td[k], 1.0}, {layout_.z[k], -instance_.endurance}},
                      GRB_LESS_EQUAL, 0.0, "endurance_" + std::to_string(k));
    }

    for (int k = 0; k < k_count; ++k) {
        for (int i = 0; i < n; ++i) {
            std::uint64_t remaining = instance_.conflict_mask[i];
            while (remaining != 0U) {
                const int j = leastBitIndex(remaining);
                remaining &= remaining - 1U;
                if (j <= i) continue;
                addConstraint({{layout_.at(layout_.y, k, i), 1.0},
                               {layout_.at(layout_.y, k, j), 1.0},
                               {layout_.z[k], -1.0}},
                              GRB_LESS_EQUAL, 0.0,
                              "conflict_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        addConstraint({{layout_.ts, instance_.vessel_speed}, {layout_.alpha[0], 1.0},
                       {layout_.at(layout_.first, 0, i), -instance_.depot_distance[i]}},
                      GRB_GREATER_EQUAL, 0.0, "start_boundary_" + std::to_string(i));
    }
    for (int k = 0; k + 1 < k_count; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                addConstraint({{layout_.beta[k], 1.0}, {layout_.tb[k], instance_.vessel_speed},
                               {layout_.alpha[k + 1], 1.0},
                               {layout_.at(layout_.last, k, i), -instance_.distance[i][j]},
                               {layout_.at(layout_.first, k + 1, j), -instance_.distance[i][j]}},
                              GRB_GREATER_EQUAL, -instance_.distance[i][j],
                              "adj_boundary_" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }
    }
}

void MasterProblem::addPhaseOneArtificialPair(int row, const std::string &name) {
    const int positive = next_column_index_++;
    int row_index = row;
    double plus = 1.0;
    checkGurobi(model_, GRBaddvar(model_, 1, &row_index, &plus, 1.0, 0.0, kInfinity,
                                  GRB_CONTINUOUS, ("art_plus_" + name).c_str()),
                 "add positive artificial");
    artificial_variables_.push_back(positive);

    const int negative = next_column_index_++;
    double minus = -1.0;
    checkGurobi(model_, GRBaddvar(model_, 1, &row_index, &minus, 1.0, 0.0, kInfinity,
                                  GRB_CONTINUOUS, ("art_minus_" + name).c_str()),
                 "add negative artificial");
    artificial_variables_.push_back(negative);
}

LinearCut MasterProblem::makeOACut(int cone_type, int slot, const Point &direction) const {
    LinearCut cut;
    cut.key = "oa:" + std::to_string(cone_type) + ":" + std::to_string(slot) + ":" +
              quantizedDirection(direction.x, direction.y);
    cut.sense = GRB_LESS_EQUAL;
    const int n = layout_.n;
    const int last_slot = layout_.slots - 1;

    switch (cone_type) {
        case 0:
            cut.coefficients = {{layout_.wx[0], direction.x}, {layout_.wy[0], direction.y},
                                {layout_.ts, -instance_.vessel_speed}};
            cut.rhs = direction.x * instance_.depot.x + direction.y * instance_.depot.y;
            break;
        case 1:
            cut.coefficients = {{layout_.qx[last_slot], direction.x}, {layout_.qy[last_slot], direction.y},
                                {layout_.tf, -instance_.vessel_speed}};
            cut.rhs = direction.x * instance_.depot.x + direction.y * instance_.depot.y;
            break;
        case 2:
            cut.coefficients = {{layout_.qx[slot], direction.x}, {layout_.qy[slot], direction.y},
                                {layout_.wx[slot], -direction.x}, {layout_.wy[slot], -direction.y},
                                {layout_.td[slot], -instance_.vessel_speed}};
            cut.rhs = 0.0;
            break;
        case 3:
            cut.coefficients = {{layout_.wx[slot + 1], direction.x}, {layout_.wy[slot + 1], direction.y},
                                {layout_.qx[slot], -direction.x}, {layout_.qy[slot], -direction.y},
                                {layout_.tb[slot], -instance_.vessel_speed}};
            cut.rhs = 0.0;
            break;
        case 4:
            cut.coefficients = {{layout_.wx[slot], direction.x}, {layout_.wy[slot], direction.y},
                                {layout_.alpha[slot], -1.0}, {layout_.z[slot], instance_.big_m_distance}};
            for (int i = 0; i < n; ++i) {
                const double coefficient = -(direction.x * instance_.platform[i].x +
                                             direction.y * instance_.platform[i].y);
                cut.coefficients.emplace_back(layout_.at(layout_.first, slot, i), coefficient);
            }
            cut.rhs = instance_.big_m_distance;
            break;
        case 5:
            cut.coefficients = {{layout_.qx[slot], -direction.x}, {layout_.qy[slot], -direction.y},
                                {layout_.beta[slot], -1.0}, {layout_.z[slot], instance_.big_m_distance}};
            for (int i = 0; i < n; ++i) {
                const double coefficient = direction.x * instance_.platform[i].x +
                                           direction.y * instance_.platform[i].y;
                cut.coefficients.emplace_back(layout_.at(layout_.last, slot, i), coefficient);
            }
            cut.rhs = instance_.big_m_distance;
            break;
        default:
            throw std::runtime_error("Unknown OA cone type");
    }
    return cut;
}

void MasterProblem::addInitialOuterApproximation() {
    const std::array<Point, 4> directions{{{1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}, {0.0, -1.0}}};
    for (const auto &u : directions) {
        addCut(makeOACut(0, 0, u));
        addCut(makeOACut(1, layout_.slots - 1, u));
        for (int k = 0; k < layout_.slots; ++k) {
            addCut(makeOACut(2, k, u));
            addCut(makeOACut(4, k, u));
            addCut(makeOACut(5, k, u));
            if (k + 1 < layout_.slots) addCut(makeOACut(3, k, u));
        }
    }
}

void MasterProblem::applyFixings(const std::vector<Fixing> &fixings) {
    for (const auto &fixing : fixings) {
        const int column = fixingColumn(fixing);
        checkGurobi(model_, GRBsetdblattrelement(model_, GRB_DBL_ATTR_LB, column, fixing.value),
                     "set fixing lower bound");
        checkGurobi(model_, GRBsetdblattrelement(model_, GRB_DBL_ATTR_UB, column, fixing.value),
                     "set fixing upper bound");
    }
}

int MasterProblem::fixingColumn(const Fixing &fixing) const {
    if (fixing.slot < 0 || fixing.slot >= layout_.slots) throw std::runtime_error("Invalid fixing slot");
    switch (fixing.kind) {
        case BranchVarKind::SortieActive:
            return layout_.z[fixing.slot];
        case BranchVarKind::PlatformAssignment:
            return layout_.at(layout_.y, fixing.slot, fixing.platform);
        case BranchVarKind::FirstPlatform:
            return layout_.at(layout_.first, fixing.slot, fixing.platform);
        case BranchVarKind::LastPlatform:
            return layout_.at(layout_.last, fixing.slot, fixing.platform);
    }
    throw std::runtime_error("Unknown fixing kind");
}

MasterSolution MasterProblem::optimize() {
    checkGurobi(model_, GRBupdatemodel(model_), "update before optimize");
    checkGurobi(model_, GRBoptimize(model_), "GRBoptimize RMP");
    int status = 0;
    checkGurobi(model_, GRBgetintattr(model_, GRB_INT_ATTR_STATUS, &status), "get RMP status");

    MasterSolution solution;
    solution.phase_one = phase_one_;
    if (status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD) {
        solution.infeasible = true;
        return solution;
    }
    if (status != GRB_OPTIMAL) {
        throw std::runtime_error("RMP did not reach optimality; status=" + std::to_string(status));
    }
    solution.optimal = true;
    checkGurobi(model_, GRBgetdblattr(model_, GRB_DBL_ATTR_OBJVAL, &solution.objective), "get RMP objective");

    int num_vars = 0;
    int num_rows = 0;
    checkGurobi(model_, GRBgetintattr(model_, GRB_INT_ATTR_NUMVARS, &num_vars), "get NumVars");
    checkGurobi(model_, GRBgetintattr(model_, GRB_INT_ATTR_NUMCONSTRS, &num_rows), "get NumConstrs");
    solution.values.resize(num_vars);
    solution.duals.resize(num_rows);
    checkGurobi(model_, GRBgetdblattrarray(model_, GRB_DBL_ATTR_X, 0, num_vars, solution.values.data()), "get RMP X");
    checkGurobi(model_, GRBgetdblattrarray(model_, GRB_DBL_ATTR_PI, 0, num_rows, solution.duals.data()), "get RMP Pi");

    solution.artificial_sum = 0.0;
    for (const int index : artificial_variables_) solution.artificial_sum += solution.values[index];
    return solution;
}

void MasterProblem::switchToPhaseTwo() {
    if (!phase_one_) return;
    for (const int artificial : artificial_variables_) {
        checkGurobi(model_, GRBsetdblattrelement(model_, GRB_DBL_ATTR_UB, artificial, 0.0),
                     "fix artificial variable");
        checkGurobi(model_, GRBsetdblattrelement(model_, GRB_DBL_ATTR_OBJ, artificial, 0.0),
                     "clear artificial objective");
    }
    checkGurobi(model_, GRBsetdblattrelement(model_, GRB_DBL_ATTR_OBJ, layout_.cmax, 1.0),
                 "set Cmax objective");
    phase_one_ = false;
    checkGurobi(model_, GRBupdatemodel(model_), "switch to phase two");
}

int MasterProblem::addPattern(int slot, int pattern_id) {
    if (slot < 0 || slot >= layout_.slots) throw std::runtime_error("Invalid slot in addPattern");
    const std::uint64_t key = slotPatternKey(slot, pattern_id);
    if (slot_pattern_present_.contains(key)) return -1;
    const RoutePattern &pattern = patterns_.get(pattern_id);

    std::vector<int> rows;
    std::vector<double> coefficients;
    rows.reserve(pattern.sequence.size() + 5);
    coefficients.reserve(pattern.sequence.size() + 5);
    rows.push_back(rows_.convexity[slot]);
    coefficients.push_back(1.0);
    for (int i = 0; i < layout_.n; ++i) {
        if (((pattern.visited_mask >> i) & 1U) == 0U) continue;
        rows.push_back(rows_.visit_link[slot * layout_.n + i]);
        coefficients.push_back(1.0);
    }
    rows.push_back(rows_.first_link[slot * layout_.n + pattern.first]);
    coefficients.push_back(1.0);
    rows.push_back(rows_.last_link[slot * layout_.n + pattern.last]);
    coefficients.push_back(1.0);
    rows.push_back(rows_.length_link[slot]);
    coefficients.push_back(-pattern.internal_distance);

    const int column = next_column_index_++;
    const std::string name = "lambda_" + std::to_string(slot) + "_" + std::to_string(pattern_id);
    checkGurobi(model_, GRBaddvar(model_, static_cast<int>(rows.size()), rows.data(), coefficients.data(),
                                  0.0, 0.0, kInfinity, GRB_CONTINUOUS, name.c_str()),
                 "add route column");
    columns_.push_back(ColumnRef{slot, pattern_id, column});
    slot_pattern_present_.insert(key);
    return column;
}

int MasterProblem::addCut(const LinearCut &cut) {
    return addConstraint(cut.coefficients, cut.sense, cut.rhs, cut.key);
}

PricingDuals MasterProblem::pricingDuals(int slot, const MasterSolution &solution) const {
    if (!solution.optimal) throw std::runtime_error("Cannot price without an optimal master solution");
    PricingDuals duals;
    duals.convexity = solution.duals.at(rows_.convexity.at(slot));
    duals.visit.resize(layout_.n);
    duals.first.resize(layout_.n);
    duals.last.resize(layout_.n);
    for (int i = 0; i < layout_.n; ++i) {
        duals.visit[i] = solution.duals.at(rows_.visit_link[slot * layout_.n + i]);
        duals.first[i] = solution.duals.at(rows_.first_link[slot * layout_.n + i]);
        duals.last[i] = solution.duals.at(rows_.last_link[slot * layout_.n + i]);
    }
    duals.length = solution.duals.at(rows_.length_link.at(slot));
    return duals;
}

PricingRestrictions MasterProblem::restrictions(int slot, const std::vector<Fixing> &fixings) const {
    PricingRestrictions result;
    int z_value = -1;
    for (const auto &fixing : fixings) {
        if (fixing.slot != slot) continue;
        const std::uint64_t bit = fixing.platform >= 0 ? (std::uint64_t{1} << fixing.platform) : 0U;
        switch (fixing.kind) {
            case BranchVarKind::SortieActive:
                if (z_value >= 0 && z_value != fixing.value) result.impossible = true;
                z_value = fixing.value;
                if (fixing.value == 0) result.skip = true;
                break;
            case BranchVarKind::PlatformAssignment:
                if (fixing.value == 1) result.required_mask |= bit;
                else result.forbidden_mask |= bit;
                break;
            case BranchVarKind::FirstPlatform: