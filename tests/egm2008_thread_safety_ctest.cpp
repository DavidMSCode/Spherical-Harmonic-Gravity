#include <array>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

#include "SHG.hpp"

namespace {

struct Case {
    double r_m;
    double phi_deg;
    double lambda_deg;
    int degree;
};

struct Reference {
    double g_mag;
    double potential;
};

constexpr double kPi = 3.14159265358979323846;
constexpr double kMachineUlps = 64.0;

double deg2rad(double deg) {
    return deg * (kPi / 180.0);
}

double norm3(const std::array<double, 3>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

bool nearly_equal_machine(double a, double b, double ulps = kMachineUlps) {
    const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
    const double tol = ulps * std::numeric_limits<double>::epsilon() * scale;
    return std::abs(a - b) <= tol;
}

bool is_finite3(const std::array<double, 3>& v) {
    return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

int run_threaded_phase(bool workspace_enabled,
                       const std::vector<Case>& cases,
                       int thread_count,
                       int iterations_per_thread) {
    SHG::set_workspace_enabled(workspace_enabled);

    std::vector<Reference> refs;
    refs.reserve(cases.size());
    for (const auto& c : cases) {
        const double phi = deg2rad(c.phi_deg);
        const double lambda = deg2rad(c.lambda_deg);
        const auto g = SHG::g_EGM2008(c.r_m, phi, lambda, c.degree);
        const double u = SHG::U_EGM2008(c.r_m, phi, lambda, c.degree);
        refs.push_back({norm3(g), u});
    }

    std::atomic<int> failures{0};

    std::vector<std::thread> threads;
    threads.reserve(thread_count);

    for (int tid = 0; tid < thread_count; ++tid) {
        threads.emplace_back([&, tid]() {
            for (int it = 0; it < iterations_per_thread; ++it) {
                for (size_t k = 0; k < cases.size(); ++k) {
                    const size_t idx = (k + static_cast<size_t>(tid + it)) % cases.size();
                    const auto& c = cases[idx];
                    const auto& ref = refs[idx];

                    const double phi = deg2rad(c.phi_deg);
                    const double lambda = deg2rad(c.lambda_deg);

                    const auto g = SHG::g_EGM2008(c.r_m, phi, lambda, c.degree);
                    const double u = SHG::U_EGM2008(c.r_m, phi, lambda, c.degree);

                    if (!is_finite3(g) || !std::isfinite(u)) {
                        failures.fetch_add(1, std::memory_order_relaxed);
                        continue;
                    }

                    const double g_mag = norm3(g);
                    if (!nearly_equal_machine(g_mag, ref.g_mag) || !nearly_equal_machine(u, ref.potential)) {
                        failures.fetch_add(1, std::memory_order_relaxed);
                    }
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    return failures.load(std::memory_order_relaxed);
}

} // namespace

int main() {
    SHG::set_coefficient_loading_verbose(false);

    const std::vector<Case> cases = {
        {6378137.0, 0.0, 0.0, 180},
        {6378237.0, 38.9072, -77.0369, 360},
        {6378137.0, 45.0, 90.0, 720},
        {6378137.0, -45.0, 0.0, 720},
        {6378137.0, 90.0, 0.0, 180},
        {6388137.0, -33.8688, 151.2093, 720},
        {6378137.0, 35.6762, 139.6503, 1080},
    };

    // Warm-up coefficient loading before concurrent execution.
    {
        const auto& c = cases.front();
        const double phi = deg2rad(c.phi_deg);
        const double lambda = deg2rad(c.lambda_deg);
        (void)SHG::g_EGM2008(c.r_m, phi, lambda, c.degree);
        (void)SHG::U_EGM2008(c.r_m, phi, lambda, c.degree);
    }

    const int thread_count = std::max(4u, std::thread::hardware_concurrency());
    const int iterations_per_thread = 40;

    const int fail_on = run_threaded_phase(true, cases, thread_count, iterations_per_thread);
    const int fail_off = run_threaded_phase(false, cases, thread_count, iterations_per_thread);

    SHG::set_workspace_enabled(true); // restore default behavior

    if (fail_on == 0 && fail_off == 0) {
        std::cout << "Thread safety checks passed (threads=" << thread_count
                  << ", iterations/thread=" << iterations_per_thread << ")." << std::endl;
        return 0;
    }

    std::cerr << "Thread safety check failed: workspace-on failures=" << fail_on
              << ", workspace-off failures=" << fail_off << std::endl;
    return 1;
}
