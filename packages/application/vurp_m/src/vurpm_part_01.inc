    std::uint64_t remaining = mask;
    while (remaining != 0U) {
        const int i = leastBitIndex(remaining);
        remaining &= remaining - 1U;
        if ((instance_.conflict_mask[i] & mask) != 0U) return false;
    }
    return true;
}

std::vector<PricingCandidate> ExactPricing::price(
    const PricingDuals &duals,
    const PricingRestrictions &restrictions,
    int max_columns) const {

    if (restrictions.skip || restrictions.impossible || max_columns <= 0) return {};
    const int n = instance_.n;
    if (n > max_exact_customers_) {
        throw std::runtime_error("Exact subset pricing supports at most " +
                                 std::to_string(max_exact_customers_) +
                                 " customers; reduce the instance or raise the configured limit with sufficient memory");
    }
    if (static_cast<int>(duals.visit.size()) != n || static_cast<int>(duals.first.size()) != n ||
        static_cast<int>(duals.last.size()) != n) {
        throw std::runtime_error("Pricing dual vector has an invalid dimension");
    }

    double length_weight = duals.length;
    if (length_weight < -1e-6) {
        throw std::runtime_error("The internal-length linking dual is negative; the master row orientation is inconsistent");
    }
    if (length_weight < 0.0) length_weight = 0.0;

    const std::uint64_t all = fullMask(n);
    const std::uint64_t allowed = all & ~restrictions.forbidden_mask;
    if ((restrictions.required_mask & ~allowed) != 0U) return {};
    if (restrictions.forced_first >= 0 && ((allowed >> restrictions.forced_first) & 1U) == 0U) return {};
    if (restrictions.forced_last >= 0 && ((allowed >> restrictions.forced_last) & 1U) == 0U) return {};

    const std::size_t state_count = std::size_t{1} << n;
    std::vector<double> visit_sum(state_count, 0.0);
    for (std::size_t mask = 1; mask < state_count; ++mask) {
        const std::uint64_t lsb = static_cast<std::uint64_t>(mask) & (~static_cast<std::uint64_t>(mask) + 1U);
        const int bit = leastBitIndex(lsb);
        visit_sum[mask] = visit_sum[mask ^ lsb] + duals.visit[bit];
    }

    std::unordered_map<std::string, PricingCandidate> best_by_sequence;
    const std::size_t local_limit = static_cast<std::size_t>(std::max(max_columns * 16, max_columns));

    for (int start = 0; start < n; ++start) {
        const std::uint64_t start_bit = std::uint64_t{1} << start;
        if ((allowed & start_bit) == 0U) continue;
        if ((restrictions.forbidden_first_mask & start_bit) != 0U) continue;
        if (restrictions.forced_first >= 0 && restrictions.forced_first != start) continue;

        const std::size_t total_entries = state_count * static_cast<std::size_t>(n);
        std::vector<double> dist(total_entries, kInfinity);
        std::vector<int> predecessor(total_entries, -2);
        const std::size_t initial_index = static_cast<std::size_t>(start_bit) * n + start;
        dist[initial_index] = 0.0;
        predecessor[initial_index] = -1;

        std::priority_queue<CandidateState, std::vector<CandidateState>, WorseCandidate> best_states;

        for (std::uint64_t mask = start_bit; mask < state_count; ++mask) {
            if ((mask & start_bit) == 0U) continue;
            if ((mask & ~allowed) != 0U) continue;
            if (!subsetConflictFree(mask)) continue;

            for (int end = 0; end < n; ++end) {
                const std::uint64_t end_bit = std::uint64_t{1} << end;
                if ((mask & end_bit) == 0U) continue;
                const std::size_t index = static_cast<std::size_t>(mask) * n + end;
                const double current_distance = dist[index];
                if (!std::isfinite(current_distance) || current_distance > instance_.max_internal_distance() + kTolerance) {
                    continue;
                }

                const bool required_satisfied = (mask & restrictions.required_mask) == restrictions.required_mask;
                const bool last_allowed = (restrictions.forbidden_last_mask & end_bit) == 0U &&
                                          (restrictions.forced_last < 0 || restrictions.forced_last == end);
                if (required_satisfied && last_allowed) {
                    const double reduced_cost = -duals.convexity - visit_sum[mask] - duals.first[start] -
                                                duals.last[end] + length_weight * current_distance;
                    if (reduced_cost < kReducedCostTolerance) {
                        CandidateState state{reduced_cost, mask, end, start, current_distance};
                        if (best_states.size() < local_limit) {
                            best_states.push(state);
                        } else if (reduced_cost < best_states.top().reduced_cost) {
                            best_states.pop();
                            best_states.push(state);
                        }
                    }
                }

                const std::uint64_t candidates = allowed & ~mask;
                std::uint64_t remaining = candidates;
                while (remaining != 0U) {
                    const int next = leastBitIndex(remaining);
                    const std::uint64_t next_bit = std::uint64_t{1} << next;
                    remaining &= remaining - 1U;
                    if ((instance_.conflict_mask[next] & mask) != 0U) continue;
                    const double next_distance = current_distance + instance_.distance[end][next];
                    if (next_distance > instance_.max_internal_distance() + kTolerance) continue;
                    const std::uint64_t next_mask = mask | next_bit;
                    const std::size_t next_index = static_cast<std::size_t>(next_mask) * n + next;
                    if (next_distance + 1e-12 < dist[next_index]) {
                        dist[next_index] = next_distance;
                        predecessor[next_index] = end;
                    }
                }
            }
        }

        std::vector<CandidateState> states;
        states.reserve(best_states.size());
        while (!best_states.empty()) {
            states.push_back(best_states.top());
            best_states.pop();
        }
        std::sort(states.begin(), states.end(), [](const CandidateState &a, const CandidateState &b) {
            return a.reduced_cost < b.reduced_cost;
        });

        for (const auto &state : states) {
            std::vector<int> reversed;
            std::uint64_t mask = state.mask;
            int end = state.end;
            while (end >= 0) {
                reversed.push_back(end);
                const std::size_t index = static_cast<std::size_t>(mask) * n + end;
                const int previous = predecessor[index];
                mask &= ~(std::uint64_t{1} << end);
                end = previous;
            }
            std::reverse(reversed.begin(), reversed.end());
            if (reversed.empty() || reversed.front() != start) continue;

            RoutePattern pattern;
            pattern.sequence = reversed;
            pattern.visited_mask = state.mask;
            pattern.first = reversed.front();
            pattern.last = reversed.back();
            pattern.internal_distance = state.distance;
            const std::string key = pattern.key();
            auto it = best_by_sequence.find(key);
            if (it == best_by_sequence.end() || state.reduced_cost < it->second.reduced_cost) {
                best_by_sequence[key] = PricingCandidate{pattern, state.reduced_cost};
            }
        }
    }

    std::vector<PricingCandidate> result;
    result.reserve(best_by_sequence.size());
    for (auto &[key, candidate] : best_by_sequence) result.push_back(std::move(candidate));
    std::sort(result.begin(), result.end(), [](const PricingCandidate &a, const PricingCandidate &b) {
        return a.reduced_cost < b.reduced_cost;
    });
    if (static_cast<int>(result.size()) > max_columns) result.resize(max_columns);
    return result;
}

std::optional<RoutePattern> ExactPricing::shortestHamiltonPath(
    std::uint64_t customer_mask,
    int first,
    int last) const {

    if (customer_mask == 0U || first < 0 || last < 0 || first >= instance_.n || last >= instance_.n) return std::nullopt;
    if ((customer_mask & (std::uint64_t{1} << first)) == 0U ||
        (customer_mask & (std::uint64_t{1} << last)) == 0U) return std::nullopt;
    if (!subsetConflictFree(customer_mask)) return std::nullopt;

    std::vector<int> nodes;
    for (int i = 0; i < instance_.n; ++i) {
        if ((customer_mask >> i) & 1U) nodes.push_back(i);
    }
    const int m = static_cast<int>(nodes.size());
    if (m > max_exact_customers_) throw std::runtime_error("Held-Karp oracle exceeds configured customer limit");
    if (m == 1) {
        if (first != last) return std::nullopt;
        RoutePattern pattern;
        pattern.sequence = {first};
        pattern.visited_mask = customer_mask;
        pattern.first = first;
        pattern.last = last;
        pattern.internal_distance = 0.0;
        return pattern;
    }
    if (first == last) return std::nullopt;

    int local_first = -1;
    int local_last = -1;
    for (int i = 0; i < m; ++i) {
        if (nodes[i] == first) local_first = i;
        if (nodes[i] == last) local_last = i;
    }
    if (local_first < 0 || local_last < 0) return std::nullopt;

    const std::size_t states = std::size_t{1} << m;
    std::vector<double> dist(states * m, kInfinity);
    std::vector<int> predecessor(states * m, -2);
    const std::uint64_t first_bit = std::uint64_t{1} << local_first;
    dist[static_cast<std::size_t>(first_bit) * m + local_first] = 0.0;
    predecessor[static_cast<std::size_t>(first_bit) * m + local_first] = -1;

    for (std::uint64_t mask = first_bit; mask < states; ++mask) {
        if ((mask & first_bit) == 0U) continue;
        for (int end = 0; end < m; ++end) {
            if (((mask >> end) & 1U) == 0U) continue;
            const double current = dist[static_cast<std::size_t>(mask) * m + end];
            if (!std::isfinite(current)) continue;
            for (int next = 0; next < m; ++next) {
                if ((mask >> next) & 1U) continue;
                if (next == local_last && popcount(mask) + 1 != m) continue;
                const std::uint64_t next_mask = mask | (std::uint64_t{1} << next);
                const double candidate = current + instance_.distance[nodes[end]][nodes[next]];
                const std::size_t next_index = static_cast<std::size_t>(next_mask) * m + next;
                if (candidate + 1e-12 < dist[next_index]) {
                    dist[next_index] = candidate;
                    predecessor[next_index] = end;
                }
            }
        }
    }

    const std::uint64_t final_mask = (std::uint64_t{1} << m) - 1U;
    const double final_distance = dist[static_cast<std::size_t>(final_mask) * m + local_last];
    if (!std::isfinite(final_distance) || final_distance > instance_.max_internal_distance() + kTolerance) {
        return std::nullopt;
    }

    std::vector<int> reversed;
    std::uint64_t mask = final_mask;
    int end = local_last;
    while (end >= 0) {
        reversed.push_back(nodes[end]);
        const int previous = predecessor[static_cast<std::size_t>(mask) * m + end];
        mask &= ~(std::uint64_t{1} << end);
        end = previous;
    }
    std::reverse(reversed.begin(), reversed.end());

    RoutePattern pattern;
    pattern.sequence = reversed;
    pattern.visited_mask = customer_mask;
    pattern.first = first;
    pattern.last = last;
    pattern.internal_distance = final_distance;
    return pattern;
}

MasterProblem::MasterProblem(
    GRBenv *environment,
    const Instance &instance,
    const PatternPool &patterns,
    const std::vector<std::vector<int>> &slot_patterns,
    const std::vector<LinearCut> &global_cuts,
    const std::vector<Fixing> &fixings,
    bool verbose_lp)
    : environment_(environment), instance_(instance), patterns_(patterns), verbose_lp_(verbose_lp) {

    if (environment_ == nullptr) throw std::runtime_error("MasterProblem received a null Gurobi environment");
    checkGurobiEnv(environment_, GRBnewmodel(environment_, &model_, "VURP_M_RMP", 0,
                                             nullptr, nullptr, nullptr, nullptr, nullptr),
                    "GRBnewmodel(VURP_M_RMP)");
    checkGurobi(model_, GRBsetintparam(GRBgetenv(model_), GRB_INT_PAR_OUTPUTFLAG, verbose_lp_ ? 1 : 0),
                 "set OutputFlag");
    checkGurobi(model_, GRBsetintparam(GRBgetenv(model_), GRB_INT_PAR_METHOD, GRB_METHOD_DUAL),
                 "set LP method");
    checkGurobi(model_, GRBsetintattr(model_, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE), "set model sense");

    layout_.n = instance_.n;
    layout_.slots = instance_.n;
    buildVariables();
    checkGurobi(model_, GRBupdatemodel(model_), "update variables");
    buildBaseRows();
    checkGurobi(model_, GRBupdatemodel(model_), "update base rows");

    for (const auto &[row, name] : phase_one_rows_) addPhaseOneArtificialPair(row, name);
    addInitialOuterApproximation();
    for (const auto &cut : global_cuts) addCut(cut);
    applyFixings(fixings);

    for (int slot = 0; slot < static_cast<int>(slot_patterns.size()); ++slot) {
        for (const int pattern_id : slot_patterns[slot]) addPattern(slot, pattern_id);
    }
    checkGurobi(model_, GRBupdatemodel(model_), "finalize RMP");
}

MasterProblem::~MasterProblem() {
    if (model_ != nullptr) GRBfreemodel(model_);
}

int MasterProblem::addVariable(double objective, double lower, double upper, char type, const std::string &name) {
    const int index = next_column_index_++;
    const std::string solver_name = gurobiName(name);
    checkGurobi(model_, GRBaddvar(model_, 0, nullptr, nullptr, objective, lower, upper, type, solver_name.c_str()),
                 "GRBaddvar(" + name + ")");
    return index;
}

int MasterProblem::addConstraint(const std::vector<std::pair<int, double>> &terms,
                                 char sense,
                                 double rhs,
                                 const std::string &name) {
    std::vector<int> indices;
    std::vector<double> coefficients;
    indices.reserve(terms.size());
    coefficients.reserve(terms.size());
    for (const auto &[index, coefficient] : terms) {
        if (std::abs(coefficient) <= 1e-14) continue;
        indices.push_back(index);
        coefficients.push_back(coefficient);
    }
    const int row = next_row_index_++;
    const std::string solver_name = gurobiName(name);
    checkGurobi(model_, GRBaddconstr(model_, static_cast<int>(indices.size()),
                                     indices.empty() ? nullptr : indices.data(),
                                     coefficients.empty() ? nullptr : coefficients.data(),
                                     sense, rhs, solver_name.c_str()),
                 "GRBaddconstr(" + name + ")");
    return row;
}

void MasterProblem::buildVariables() {
    const int n = layout_.n;
    const int k_count = layout_.slots;
    layout_.z.resize(k_count);
    layout_.ell.resize(k_count);
    layout_.wx.resize(k_count);
    layout_.wy.resize(k_count);
    layout_.qx.resize(k_count);
    layout_.qy.resize(k_count);
    layout_.alpha.resize(k_count);
    layout_.beta.resize(k_count);
    layout_.td.resize(k_count);
    layout_.tb.resize(std::max(0, k_count - 1));
    layout_.y.resize(k_count * n);
    layout_.first.resize(k_count * n);
    layout_.last.resize(k_count * n);

    for (int k = 0; k < k_count; ++k) {
        layout_.z[k] = addVariable(0.0, 0.0, 1.0, GRB_CONTINUOUS, "z_" + std::to_string(k));
    }
    for (int k = 0; k < k_count; ++k) {
        for (int i = 0; i < n; ++i) {
            layout_.y[k * n + i] = addVariable(0.0, 0.0, 1.0, GRB_CONTINUOUS,
                                               "y_" + std::to_string(i) + "_" + std::to_string(k));
            layout_.first[k * n + i] = addVariable(0.0, 0.0, 1.0, GRB_CONTINUOUS,
                                                   "first_" + std::to_string(i) + "_" + std::to_string(k));
            layout_.last[k * n + i] = addVariable(0.0, 0.0, 1.0, GRB_CONTINUOUS,
                                                  "last_" + std::to_string(i) + "_" + std::to_string(k));
        }
    }
    for (int k = 0; k < k_count; ++k) {
        layout_.ell[k] = addVariable(0.0, 0.0, instance_.max_internal_distance(), GRB_CONTINUOUS,
                                     "ell_" + std::to_string(k));
        layout_.wx[k] = addVariable(0.0, instance_.omega.xmin, instance_.omega.xmax, GRB_CONTINUOUS,
                                    "wx_" + std::to_string(k));
        layout_.wy[k] = addVariable(0.0, instance_.omega.ymin, instance_.omega.ymax, GRB_CONTINUOUS,
