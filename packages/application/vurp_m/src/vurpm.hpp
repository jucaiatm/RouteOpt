#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gurobi_c.h"

namespace RouteOpt::Application::VURPM {

constexpr double kTolerance = 1e-7;
constexpr double kReducedCostTolerance = -1e-7;
constexpr double kIntegralityTolerance = 1e-6;
constexpr double kOATolerance = 1e-7;
constexpr double kInfinity = GRB_INFINITY;

struct Point {
    double x{};
    double y{};
};

struct OmegaBox {
    double xmin{};
    double xmax{};
    double ymin{};
    double ymax{};

    [[nodiscard]] bool contains(const Point &p, double tol = kTolerance) const;
};

struct Instance {
    std::string name{"VURP-M"};
    int n{};
    double vessel_speed{};
    double uav_speed{};
    double endurance{};
    Point depot{};
    OmegaBox omega{};
    bool omega_explicit{};
    std::vector<Point> platform{};
    std::vector<std::vector<double>> distance{};
    std::vector<double> depot_distance{};
    std::vector<std::uint64_t> conflict_mask{};
    double big_m_distance{};

    static Instance read(const std::string &path);
    void preprocess();
    void validate() const;

    [[nodiscard]] double max_internal_distance() const {
        return uav_speed * endurance;
    }
};

struct RoutePattern {
    int id{-1};
    std::vector<int> sequence{};
    std::uint64_t visited_mask{};
    int first{-1};
    int last{-1};
    double internal_distance{};

    [[nodiscard]] std::string key() const;
};

class PatternPool {
public:
    int add(const RoutePattern &pattern);
    [[nodiscard]] const RoutePattern &get(int id) const;
    [[nodiscard]] int size() const;
    [[nodiscard]] std::optional<int> find(const RoutePattern &pattern) const;

private:
    std::vector<RoutePattern> patterns_{};
    std::unordered_map<std::string, int> index_{};
};

enum class BranchVarKind {
    SortieActive,
    PlatformAssignment,
    FirstPlatform,
    LastPlatform
};

struct Fixing {
    BranchVarKind kind{BranchVarKind::PlatformAssignment};
    int slot{-1};
    int platform{-1};
    int value{};
};

struct PricingRestrictions {
    bool skip{};
    bool impossible{};
    std::uint64_t required_mask{};
    std::uint64_t forbidden_mask{};
    std::uint64_t forbidden_first_mask{};
    std::uint64_t forbidden_last_mask{};
    int forced_first{-1};
    int forced_last{-1};
};

struct PricingDuals {
    double convexity{};
    std::vector<double> visit{};
    std::vector<double> first{};
    std::vector<double> last{};
    double length{};
};

struct PricingCandidate {
    RoutePattern pattern{};
    double reduced_cost{};
};

class ExactPricing {
public:
    explicit ExactPricing(const Instance &instance, int max_exact_customers = 20);

    [[nodiscard]] std::vector<PricingCandidate> price(
        const PricingDuals &duals,
        const PricingRestrictions &restrictions,
        int max_columns) const;

    [[nodiscard]] std::optional<RoutePattern> shortestHamiltonPath(
        std::uint64_t customer_mask,
        int first,
        int last) const;

private:
    const Instance &instance_;
    int max_exact_customers_{};

    [[nodiscard]] bool subsetConflictFree(std::uint64_t mask) const;
};

struct LinearCut {
    std::string key{};
    std::vector<std::pair<int, double>> coefficients{};
    char sense{GRB_LESS_EQUAL};
    double rhs{};
};

struct MasterLayout {
    int n{};
    int slots{};
    std::vector<int> z{};
    std::vector<int> y{};
    std::vector<int> first{};
    std::vector<int> last{};
    std::vector<int> ell{};
    std::vector<int> wx{};
    std::vector<int> wy{};
    std::vector<int> qx{};
    std::vector<int> qy{};
    std::vector<int> alpha{};
    std::vector<int> beta{};
    std::vector<int> td{};
    std::vector<int> tb{};
    int ts{-1};
    int tf{-1};
    int cmax{-1};
    int static_variable_count{};

    [[nodiscard]] int at(const std::vector<int> &v, int slot, int platform) const {
        return v[slot * n + platform];
    }
};

struct MasterRows {
    std::vector<int> convexity{};
    std::vector<int> visit_link{};
    std::vector<int> first_link{};
    std::vector<int> last_link{};
    std::vector<int> length_link{};
};

struct ColumnRef {
    int slot{-1};
    int pattern_id{-1};
    int column_index{-1};
};

struct MasterSolution {
    bool optimal{};
    bool infeasible{};
    bool phase_one{};
    double objective{kInfinity};
    double artificial_sum{kInfinity};
    std::vector<double> values{};
    std::vector<double> duals{};
};

struct BranchChoice {
    Fixing zero_child{};
    Fixing one_child{};
    double fractionality{};
};

class MasterProblem {
public:
    MasterProblem(
        GRBenv *environment,
        const Instance &instance,
        const PatternPool &patterns,
        const std::vector<std::vector<int>> &slot_patterns,
        const std::vector<LinearCut> &global_cuts,
        const std::vector<Fixing> &fixings,
        bool verbose_lp = false);

    ~MasterProblem();
    MasterProblem(const MasterProblem &) = delete;
    MasterProblem &operator=(const MasterProblem &) = delete;

    MasterSolution optimize();
    void switchToPhaseTwo();
    int addPattern(int slot, int pattern_id);
    int addCut(const LinearCut &cut);

    [[nodiscard]] PricingDuals pricingDuals(int slot, const MasterSolution &solution) const;
    [[nodiscard]] PricingRestrictions restrictions(int slot, const std::vector<Fixing> &fixings) const;
    [[nodiscard]] std::vector<LinearCut> separateOuterApproximation(
        const MasterSolution &solution,
        int maximum_cuts) const;
    [[nodiscard]] bool staticVariablesIntegral(const MasterSolution &solution) const;
    [[nodiscard]] std::optional<BranchChoice> chooseBranch(const MasterSolution &solution) const;
    [[nodiscard]] std::string partitionKey(const MasterSolution &solution) const;
    [[nodiscard]] const MasterLayout &layout() const { return layout_; }
    [[nodiscard]] GRBmodel *rawModel() const { return model_; }
    [[nodiscard]] bool inPhaseOne() const { return phase_one_; }

private:
    GRBenv *environment_{};
    GRBmodel *model_{};
    const Instance &instance_;
    const PatternPool &patterns_;
    MasterLayout layout_{};
    MasterRows rows_{};
    bool phase_one_{true};
    bool verbose_lp_{};
    std::vector<int> artificial_variables_{};
    std::vector<ColumnRef> columns_{};
    std::unordered_set<std::uint64_t> slot_pattern_present_{};
    std::vector<std::pair<int, std::string>> phase_one_rows_{};
    int next_column_index_{};
    int next_row_index_{};

    int addVariable(double objective, double lower, double upper, char type, const std::string &name);
    int addConstraint(const std::vector<std::pair<int, double>> &terms, char sense, double rhs,
                      const std::string &name);
    void buildVariables();
    void buildBaseRows();
    void addInitialOuterApproximation();
    void addPhaseOneArtificialPair(int row, const std::string &name);
    void applyFixings(const std::vector<Fixing> &fixings);
    [[nodiscard]] int fixingColumn(const Fixing &fixing) const;
    [[nodiscard]] LinearCut makeOACut(int cone_type, int slot, const Point &direction) const;
};

struct SortiePlan {
    int slot{-1};
    std::uint64_t customer_mask{};
    int first{-1};
    int last{-1};
    RoutePattern route{};
};

struct PartitionPlan {
    std::vector<SortiePlan> sorties{};
    std::vector<int> z{};
    std::vector<int> y{};
    std::vector<int> first{};
    std::vector<int> last{};
    std::string key{};
};

struct SocpSolution {
    bool feasible{};
    bool proven_infeasible{};
    double objective{kInfinity};
    std::vector<Point> launch{};
    std::vector<Point> recovery{};
    std::vector<double> sortie_time{};
    std::vector<double> bridge_time{};
    double start_time{};
    double finish_time{};
};

class SocpOracle {
public:
    SocpOracle(GRBenv *environment, const Instance &instance, bool verbose = false);
    [[nodiscard]] SocpSolution solve(const PartitionPlan &partition) const;

private:
    GRBenv *environment_{};
    const Instance &instance_;
    bool verbose_{};
};

struct BpcParameters {
    double time_limit{7200.0};
    int max_pricing_columns_per_slot{20};
    int max_oa_cuts_per_round{100};
    int max_exact_pricing_customers{20};
    int threads{1};
    bool verbose_lp{};
    bool verbose_socp{};
};

struct Incumbent {
    bool available{};
    double objective{kInfinity};
    PartitionPlan partition{};
    SocpSolution continuous{};
};

struct BpcStatistics {
    int explored_nodes{};
    int generated_patterns{};
    int generated_slot_columns{};
    int oa_cuts{};
    int logic_cuts{};
    int socp_calls{};
    double lower_bound{};
    double elapsed{};
};

class BpcSolver {
public:
    BpcSolver(Instance instance, BpcParameters parameters);
    ~BpcSolver();
    BpcSolver(const BpcSolver &) = delete;
    BpcSolver &operator=(const BpcSolver &) = delete;

    void solve();
    void printResult() const;

private:
    struct Node {
        int id{};
        double estimate{};
        std::vector<Fixing> fixings{};
    };

    struct NodeCompare {
        bool operator()(const Node &a, const Node &b) const {
            if (std::abs(a.estimate - b.estimate) > kTolerance) return a.estimate > b.estimate;
            return a.id > b.id;
        }
    };

    struct PartitionEvaluation {
        bool feasible{};
        bool proven_infeasible{};
        double value{kInfinity};
        PartitionPlan partition{};
        SocpSolution socp{};
    };

    Instance instance_;
    BpcParameters parameters_;
    GRBenv *environment_{};
    PatternPool patterns_{};
    std::vector<std::vector<int>> slot_patterns_{};
    std::vector<LinearCut> global_cuts_{};
    std::unordered_set<std::string> global_cut_keys_{};
    std::unordered_map<std::string, PartitionEvaluation> partition_cache_{};
    ExactPricing pricing_;
    std::unique_ptr<SocpOracle> socp_oracle_{};
    Incumbent incumbent_{};
    BpcStatistics statistics_{};
    std::chrono::steady_clock::time_point start_time_{};
    int next_node_id_{};

    void initializeEnvironment();
    void initializePatterns();
    void initializeIncumbent();
    [[nodiscard]] bool timeLimitReached() const;
    [[nodiscard]] double elapsed() const;
    [[nodiscard]] bool fixingSetContradictory(const std::vector<Fixing> &fixings) const;
    [[nodiscard]] PartitionPlan extractPartition(const MasterProblem &master,
                                                 const MasterSolution &solution) const;
    [[nodiscard]] LinearCut buildPartitionValueCut(const MasterLayout &layout,
                                                   const PartitionEvaluation &evaluation) const;
    [[nodiscard]] LinearCut buildPartitionNoGoodCut(const MasterLayout &layout,
                                                    const PartitionPlan &partition) const;
    bool registerGlobalCut(const LinearCut &cut);
    int registerPatternForSlot(int slot, RoutePattern pattern, MasterProblem &master);
};

} // namespace RouteOpt::Application::VURPM
