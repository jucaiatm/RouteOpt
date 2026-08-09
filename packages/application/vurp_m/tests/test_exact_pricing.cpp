#include "vurpm.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using RouteOpt::Application::VURPM::ExactPricing;
using RouteOpt::Application::VURPM::Instance;
using RouteOpt::Application::VURPM::OmegaBox;
using RouteOpt::Application::VURPM::Point;
using RouteOpt::Application::VURPM::PricingDuals;
using RouteOpt::Application::VURPM::PricingRestrictions;
using RouteOpt::Application::VURPM::RoutePattern;

namespace {

constexpr double kTestTolerance = 1e-9;

struct BruteResult {
    bool found{};
    double reduced_cost{std::numeric_limits<double>::infinity()};
    std::vector<int> sequence{};
    double internal_distance{};
};

std::uint64_t allMask(int n) {
    return (std::uint64_t{1} << n) - 1U;
}

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

bool conflictFree(const Instance &instance, std::uint64_t mask) {
    for (int i = 0; i < instance.n; ++i) {
        if (((mask >> i) & 1U) == 0U) continue;
        if ((instance.conflict_mask[i] & mask) != 0U) return false;
    }
    return true;
}

BruteResult brutePrice(const Instance &instance,
                       const PricingDuals &duals,
                       const PricingRestrictions &restrictions) {
    BruteResult best;
    if (restrictions.skip || restrictions.impossible) return best;

    const std::uint64_t allowed = allMask(instance.n) & ~restrictions.forbidden_mask;
    if ((restrictions.required_mask & ~allowed) != 0U) return best;

    std::vector<int> path;
    std::function<void(int, int, std::uint64_t, double, double)> dfs;
    dfs = [&](int start, int end, std::uint64_t mask, double distance, double visit_sum) {
        const std::uint64_t end_bit = std::uint64_t{1} << end;
        const bool required_satisfied =
            (mask & restrictions.required_mask) == restrictions.required_mask;
        const bool last_allowed =
            (restrictions.forbidden_last_mask & end_bit) == 0U &&
            (restrictions.forced_last < 0 || restrictions.forced_last == end);

        if (required_satisfied && last_allowed) {
            const double reduced_cost =
                -duals.convexity - visit_sum - duals.first[start] - duals.last[end] +
                duals.length * distance;
            if (!best.found || reduced_cost < best.reduced_cost - kTestTolerance) {
                best.found = true;
                best.reduced_cost = reduced_cost;
                best.sequence = path;
                best.internal_distance = distance;
            }
        }

        for (int next = 0; next < instance.n; ++next) {
            const std::uint64_t next_bit = std::uint64_t{1} << next;
            if ((allowed & next_bit) == 0U || (mask & next_bit) != 0U) continue;
            if ((instance.conflict_mask[next] & mask) != 0U) continue;
            const double next_distance = distance + instance.distance[end][next];
            if (next_distance > instance.max_internal_distance() + kTestTolerance) continue;
            path.push_back(next);
            dfs(start, next, mask | next_bit, next_distance,
                visit_sum + duals.visit[next]);
            path.pop_back();
        }
    };

    for (int start = 0; start < instance.n; ++start) {
        const std::uint64_t bit = std::uint64_t{1} << start;
        if ((allowed & bit) == 0U) continue;
        if ((restrictions.forbidden_first_mask & bit) != 0U) continue;
        if (restrictions.forced_first >= 0 && restrictions.forced_first != start) continue;
        path.assign(1, start);
        dfs(start, start, bit, 0.0, duals.visit[start]);
    }
    return best;
}

double bruteHamiltonDistance(const Instance &instance,
                            std::uint64_t mask,
                            int first,
                            int last) {
    std::vector<int> middle;
    for (int i = 0; i < instance.n; ++i) {
        if (((mask >> i) & 1U) != 0U && i != first && i != last) middle.push_back(i);
    }
    std::sort(middle.begin(), middle.end());
    double best = std::numeric_limits<double>::infinity();
    do {
        double distance = 0.0;
        int previous = first;
        for (const int node : middle) {
            distance += instance.distance[previous][node];
            previous = node;
        }
        distance += instance.distance[previous][last];
        best = std::min(best, distance);
    } while (std::next_permutation(middle.begin(), middle.end()));
    return best;
}

Instance makeInstance() {
    Instance instance;
    instance.name = "pricing-regression";
    instance.n = 5;
    instance.vessel_speed = 1.0;
    instance.uav_speed = 3.0;
    instance.endurance = 10.0;
    instance.depot = Point{0.0, 0.0};
    instance.omega = OmegaBox{-2.0, 5.0, -2.0, 5.0};
    instance.omega_explicit = true;
    instance.platform = {
        Point{1.0, 0.0},
        Point{2.0, 0.0},
        Point{1.0, 1.0},
        Point{2.0, 1.0},
        Point{1.5, 0.5},
    };
    instance.preprocess();
    instance.validate();
    return instance;
}

void comparePricing(const Instance &instance,
                    const ExactPricing &pricing,
                    const PricingDuals &duals,
                    const PricingRestrictions &restrictions,
                    const std::string &case_name) {
    const BruteResult brute = brutePrice(instance, duals, restrictions);
    const auto generated = pricing.price(duals, restrictions, 512);

    if (!brute.found || brute.reduced_cost >= -1e-7) {
        require(generated.empty(), case_name + ": pricing returned a false negative-cost column");
        return;
    }

    require(!generated.empty(), case_name + ": exact pricing missed a negative-cost column");
    require(std::abs(generated.front().reduced_cost - brute.reduced_cost) <= 1e-8,
            case_name + ": exact pricing reduced cost differs from brute force");
    require(conflictFree(instance, generated.front().pattern.visited_mask),
            case_name + ": pricing generated an incompatible platform subset");
}

} // namespace

int main() {
    try {
        const Instance instance = makeInstance();
        const ExactPricing pricing(instance, 20);

        PricingDuals duals;
        duals.convexity = 0.2;
        duals.visit = {1.30, 0.80, 1.10, 0.70, 0.90};
        duals.first = {0.40, 0.30, 0.20, 0.10, 0.50};
        duals.last = {0.20, 0.40, 0.30, 0.50, 0.10};
        duals.length = 0.15;

        comparePricing(instance, pricing, duals, PricingRestrictions{}, "unrestricted");

        PricingRestrictions restricted;
        restricted.required_mask = std::uint64_t{1} << 2;
        restricted.forbidden_mask = std::uint64_t{1} << 4;
        restricted.forced_first = 1;
        restricted.forced_last = 3;
        comparePricing(instance, pricing, duals, restricted, "branch-restricted");

        PricingRestrictions impossible = restricted;
        impossible.forbidden_mask |= std::uint64_t{1} << 2;
        require(pricing.price(duals, impossible, 32).empty(),
                "inconsistent required/forbidden restrictions were not rejected");

        PricingDuals zero_duals;
        zero_duals.visit.assign(instance.n, 0.0);
        zero_duals.first.assign(instance.n, 0.0);
        zero_duals.last.assign(instance.n, 0.0);
        zero_duals.length = 0.0;
        require(pricing.price(zero_duals, PricingRestrictions{}, 32).empty(),
                "zero-dual pricing should not generate a negative reduced-cost column");

        const std::uint64_t hamilton_mask =
            (std::uint64_t{1} << 0) |
            (std::uint64_t{1} << 1) |
            (std::uint64_t{1} << 2) |
            (std::uint64_t{1} << 3);
        const auto route = pricing.shortestHamiltonPath(hamilton_mask, 0, 3);
        require(route.has_value(), "Held--Karp oracle failed on a feasible fixed-endpoint path");
        const double brute_distance =
            bruteHamiltonDistance(instance, hamilton_mask, 0, 3);
        require(std::abs(route->internal_distance - brute_distance) <= 1e-9,
                "Held--Karp distance differs from permutation enumeration");
        require(route->sequence.front() == 0 && route->sequence.back() == 3,
                "Held--Karp route violates fixed endpoints");

        std::cout << "VURP-M exact-pricing regression passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "VURP-M exact-pricing regression failed: " << error.what() << '\n';
        return 1;
    }
}
