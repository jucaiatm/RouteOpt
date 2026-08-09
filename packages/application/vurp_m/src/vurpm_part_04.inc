        }
        if (fixing.kind == BranchVarKind::LastPlatform && fixing.value == 1) {
            if (forced_last[fixing.slot] >= 0 && forced_last[fixing.slot] != fixing.platform) return true;
            forced_last[fixing.slot] = fixing.platform;
            required[fixing.slot] |= std::uint64_t{1} << fixing.platform;
        }
    }
    for (int k = 0; k < instance_.n; ++k) {
        if ((required[k] & forbidden[k]) != 0U) return true;
        if (fixed_z[k] == 0 && required[k] != 0U) return true;
        std::uint64_t r = required[k];
        while (r != 0U) {
            const int i = leastBitIndex(r);
            r &= r - 1U;
            if ((instance_.conflict_mask[i] & required[k]) != 0U) return true;
        }
    }
    return false;
}

PartitionPlan BpcSolver::extractPartition(const MasterProblem &master,
                                          const MasterSolution &solution) const {
    const MasterLayout &layout = master.layout();
    PartitionPlan partition;
    partition.key = master.partitionKey(solution);
    partition.z.resize(layout.slots);
    partition.y.resize(layout.slots * layout.n);
    partition.first.resize(layout.slots * layout.n);
    partition.last.resize(layout.slots * layout.n);

    for (int k = 0; k < layout.slots; ++k) partition.z[k] = solution.values[layout.z[k]] >= 0.5 ? 1 : 0;
    for (int pos = 0; pos < layout.slots * layout.n; ++pos) {
        partition.y[pos] = solution.values[layout.y[pos]] >= 0.5 ? 1 : 0;
        partition.first[pos] = solution.values[layout.first[pos]] >= 0.5 ? 1 : 0;
        partition.last[pos] = solution.values[layout.last[pos]] >= 0.5 ? 1 : 0;
    }

    for (int k = 0; k < layout.slots; ++k) {
        if (partition.z[k] == 0) continue;
        std::uint64_t mask = 0U;
        int first = -1;
        int last = -1;
        for (int i = 0; i < layout.n; ++i) {
            if (partition.y[k * layout.n + i]) mask |= std::uint64_t{1} << i;
            if (partition.first[k * layout.n + i]) first = i;
            if (partition.last[k * layout.n + i]) last = i;
        }
        const auto route = pricing_.shortestHamiltonPath(mask, first, last);
        if (!route) {
            SortiePlan invalid;
            invalid.slot = k;
            invalid.customer_mask = mask;
            invalid.first = first;
            invalid.last = last;
            partition.sorties.push_back(invalid);
        } else {
            partition.sorties.push_back(SortiePlan{k, mask, first, last, *route});
        }
    }
    return partition;
}

LinearCut BpcSolver::buildPartitionValueCut(const MasterLayout &layout,
                                            const PartitionEvaluation &evaluation) const {
    const double q = evaluation.value;
    LinearCut cut;
    cut.key = "partition_value:" + evaluation.partition.key;
    cut.sense = GRB_GREATER_EQUAL;
    cut.coefficients.emplace_back(layout.cmax, 1.0);
    int ones = 0;
    auto append = [&](const std::vector<int> &indices, const std::vector<int> &values) {
        for (std::size_t i = 0; i < indices.size(); ++i) {
            if (values[i] == 1) {
                cut.coefficients.emplace_back(indices[i], -q);
                ++ones;
            } else {
                cut.coefficients.emplace_back(indices[i], q);
            }
        }
    };
    append(layout.z, evaluation.partition.z);
    append(layout.y, evaluation.partition.y);
    append(layout.first, evaluation.partition.first);
    append(layout.last, evaluation.partition.last);
    cut.rhs = q * (1.0 - static_cast<double>(ones));
    return cut;
}

LinearCut BpcSolver::buildPartitionNoGoodCut(const MasterLayout &layout,
                                             const PartitionPlan &partition) const {
    LinearCut cut;
    cut.key = "partition_nogood:" + partition.key;
    cut.sense = GRB_GREATER_EQUAL;
    int ones = 0;
    auto append = [&](const std::vector<int> &indices, const std::vector<int> &values) {
        for (std::size_t i = 0; i < indices.size(); ++i) {
            if (values[i] == 1) {
                cut.coefficients.emplace_back(indices[i], -1.0);
                ++ones;
            } else {
                cut.coefficients.emplace_back(indices[i], 1.0);
            }
        }
    };
    append(layout.z, partition.z);
    append(layout.y, partition.y);
    append(layout.first, partition.first);
    append(layout.last, partition.last);
    cut.rhs = 1.0 - static_cast<double>(ones);
    return cut;
}

bool BpcSolver::registerGlobalCut(const LinearCut &cut) {
    if (!global_cut_keys_.insert(cut.key).second) return false;
    global_cuts_.push_back(cut);
    return true;
}

int BpcSolver::registerPatternForSlot(int slot, RoutePattern pattern, MasterProblem &master) {
    const int old_size = patterns_.size();
    const int id = patterns_.add(pattern);
    if (patterns_.size() > old_size) {
        ++statistics_.generated_patterns;
    }
    auto &slot_patterns = slot_patterns_[slot];
    if (std::find(slot_patterns.begin(), slot_patterns.end(), id) != slot_patterns.end()) return 0;
    slot_patterns.push_back(id);
    master.addPattern(slot, id);
    ++statistics_.generated_slot_columns;
    return 1;
}

void BpcSolver::solve() {
    start_time_ = std::chrono::steady_clock::now();
    initializePatterns();
    initializeIncumbent();

    std::priority_queue<Node, std::vector<Node>, NodeCompare> queue;
    queue.push(Node{next_node_id_++, 0.0, {}});
    double best_open_bound = 0.0;

    while (!queue.empty() && !timeLimitReached()) {
        Node node = queue.top();
        queue.pop();
        if (incumbent_.available && node.estimate >= incumbent_.objective - kTolerance) continue;
        if (fixingSetContradictory(node.fixings)) continue;
        ++statistics_.explored_nodes;

        MasterProblem master(environment_, instance_, patterns_, slot_patterns_, global_cuts_, node.fixings,
                             parameters_.verbose_lp);
        bool node_finished = false;
        while (!node_finished && !timeLimitReached()) {
            MasterSolution solution = master.optimize();
            if (solution.infeasible) break;

            int added_columns = 0;
            for (int k = 0; k < instance_.n; ++k) {
                const PricingRestrictions restrictions = master.restrictions(k, node.fixings);
                if (restrictions.skip || restrictions.impossible) continue;
                const PricingDuals duals = master.pricingDuals(k, solution);
                const auto candidates = pricing_.price(duals, restrictions,
                                                       parameters_.max_pricing_columns_per_slot);
                for (const auto &candidate : candidates) {
                    if (candidate.reduced_cost >= kReducedCostTolerance) continue;
                    added_columns += registerPatternForSlot(k, candidate.pattern, master);
                }
            }
            if (added_columns > 0) continue;

            if (master.inPhaseOne()) {
                if (solution.artificial_sum > 1e-6) {
                    node_finished = true;
                    break;
                }
                master.switchToPhaseTwo();
                continue;
            }

            if (incumbent_.available && solution.objective >= incumbent_.objective - kTolerance) {
                node_finished = true;
                break;
            }

            int added_oa = 0;
            for (const auto &cut : master.separateOuterApproximation(solution,
                                                                     parameters_.max_oa_cuts_per_round)) {
                if (!registerGlobalCut(cut)) continue;
                master.addCut(cut);
                ++statistics_.oa_cuts;
                ++added_oa;
            }
            if (added_oa > 0) continue;

            node.estimate = solution.objective;
            best_open_bound = solution.objective;

            if (master.staticVariablesIntegral(solution)) {
                const std::string key = master.partitionKey(solution);
                auto cache_it = partition_cache_.find(key);
                if (cache_it == partition_cache_.end()) {
                    PartitionPlan partition = extractPartition(master, solution);
                    PartitionEvaluation evaluation;
                    evaluation.partition = partition;
                    bool route_valid = !partition.sorties.empty();
                    for (const auto &sortie : partition.sorties) {
                        if (sortie.route.sequence.empty()) {
                            route_valid = false;
                            break;
                        }
                    }
                    if (!route_valid) {
                        evaluation.proven_infeasible = true;
                    } else {
                        ++statistics_.socp_calls;
                        evaluation.socp = socp_oracle_->solve(partition);
                        evaluation.feasible = evaluation.socp.feasible;
                        evaluation.proven_infeasible = evaluation.socp.proven_infeasible;
                        evaluation.value = evaluation.socp.objective;
                    }
                    auto [it, inserted] = partition_cache_.emplace(key, std::move(evaluation));
                    cache_it = it;

                    LinearCut cut;
                    if (cache_it->second.feasible) {
                        if (!incumbent_.available || cache_it->second.value < incumbent_.objective - kTolerance) {
                            incumbent_.available = true;
                            incumbent_.objective = cache_it->second.value;
                            incumbent_.partition = cache_it->second.partition;
                            incumbent_.continuous = cache_it->second.socp;
                            std::cout << "Certified UB: " << std::setprecision(12) << incumbent_.objective
                                      << " at node " << node.id << '\n';
                        }
                        cut = buildPartitionValueCut(master.layout(), cache_it->second);
                    } else if (cache_it->second.proven_infeasible) {
                        cut = buildPartitionNoGoodCut(master.layout(), cache_it->second.partition);
                    } else {
                        throw std::runtime_error("SOCP oracle returned an undecided status");
                    }
                    if (registerGlobalCut(cut)) {
                        master.addCut(cut);
                        ++statistics_.logic_cuts;
                        continue;
                    }
                }

                const PartitionEvaluation &evaluation = cache_it->second;
                if (evaluation.feasible &&
                    (!incumbent_.available || evaluation.value < incumbent_.objective - kTolerance)) {
                    incumbent_.available = true;
                    incumbent_.objective = evaluation.value;
                    incumbent_.partition = evaluation.partition;
                    incumbent_.continuous = evaluation.socp;
                }
                node_finished = true;
                break;
            }

            const auto branch = master.chooseBranch(solution);
            if (!branch) {
                throw std::runtime_error("Fractional master solution has no branching candidate");
            }
            Node zero_child{next_node_id_++, solution.objective, node.fixings};
            zero_child.fixings.push_back(branch->zero_child);
            Node one_child{next_node_id_++, solution.objective, node.fixings};
            one_child.fixings.push_back(branch->one_child);
            if (!fixingSetContradictory(zero_child.fixings)) queue.push(std::move(zero_child));
            if (!fixingSetContradictory(one_child.fixings)) queue.push(std::move(one_child));
            node_finished = true;
        }
    }

    statistics_.elapsed = elapsed();
    if (queue.empty() && incumbent_.available) {
        statistics_.lower_bound = incumbent_.objective;
    } else if (!queue.empty()) {
        statistics_.lower_bound = queue.top().estimate;
    } else {
        statistics_.lower_bound = best_open_bound;
    }
}

void BpcSolver::printResult() const {
    std::cout << "\n===== VURP-M BPC result =====\n";
    std::cout << "Instance: " << instance_.name << '\n';
    std::cout << "Status: " << (incumbent_.available ? "feasible" : "no certified solution") << '\n';
    if (incumbent_.available) {
        std::cout << std::setprecision(12);
        std::cout << "Objective: " << incumbent_.objective << '\n';
        std::cout << "Lower bound: " << statistics_.lower_bound << '\n';
        const double gap = std::max(0.0, incumbent_.objective - statistics_.lower_bound) /
                           std::max(1.0, std::abs(incumbent_.objective));
        std::cout << "Relative gap: " << gap << '\n';
        std::cout << "Sorties: " << incumbent_.partition.sorties.size() << '\n';
        for (std::size_t k = 0; k < incumbent_.partition.sorties.size(); ++k) {
            const auto &sortie = incumbent_.partition.sorties[k];
            std::cout << "  " << (k + 1) << ": ";
            for (std::size_t h = 0; h < sortie.route.sequence.size(); ++h) {
                if (h != 0) std::cout << " -> ";
                std::cout << sortie.route.sequence[h] + 1;
            }
            std::cout << " | internal=" << sortie.route.internal_distance;
            if (k < incumbent_.continuous.launch.size()) {
                std::cout << " | launch=(" << incumbent_.continuous.launch[k].x << ','
                          << incumbent_.continuous.launch[k].y << ")";
                std::cout << " recovery=(" << incumbent_.continuous.recovery[k].x << ','
                          << incumbent_.continuous.recovery[k].y << ")";
                std::cout << " duration=" << incumbent_.continuous.sortie_time[k];
            }
            std::cout << '\n';
        }
    }
    std::cout << "Nodes: " << statistics_.explored_nodes << '\n';
    std::cout << "Patterns: " << statistics_.generated_patterns << '\n';
    std::cout << "Slot columns: " << statistics_.generated_slot_columns << '\n';
    std::cout << "OA cuts: " << statistics_.oa_cuts << '\n';
    std::cout << "Logic cuts: " << statistics_.logic_cuts << '\n';
    std::cout << "SOCP calls: " << statistics_.socp_calls << '\n';
    std::cout << "Elapsed: " << statistics_.elapsed << " s\n";
}

} // namespace RouteOpt::Application::VURPM
