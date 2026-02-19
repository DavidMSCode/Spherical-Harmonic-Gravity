#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "SHG.hpp"

namespace {
struct Case {
    double r_m;
    double phi_deg;
    double lambda_deg;
    int degree;
};

constexpr double kPi = 3.14159265358979323846;
constexpr double kGTol = 1e-16;
constexpr double kUTol = 1e-16;

double deg2rad(double deg) {
    return deg * (kPi / 180.0);
}

double norm3(const std::array<double, 3>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

bool nearly_equal(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

} // namespace

int main() {
    const std::vector<Case> cases = {
        {6378137.0, 0.0, 0.0, 180},
        {6378237.0, 38.9072, -77.0369, 360},
        {6378137.0, 35.6762, 139.6503, 720},
        {6388137.0, -33.8688, 151.2093, 720},
        {6378137.0, 90.0, 0.0, 180},
    };

    SHG::set_coefficient_loading_verbose(false);

    int failures = 0;

    // Numerical equivalence test: workspace ON vs OFF
    for (const auto& tc : cases) {
        const double phi = deg2rad(tc.phi_deg);
        const double lambda = deg2rad(tc.lambda_deg);

        SHG::set_workspace_enabled(true);
        const auto g_on = SHG::g_EGM2008(tc.r_m, phi, lambda, tc.degree);
        const double u_on = SHG::U_EGM2008(tc.r_m, phi, lambda, tc.degree);

        SHG::set_workspace_enabled(false);
        const auto g_off = SHG::g_EGM2008(tc.r_m, phi, lambda, tc.degree);
        const double u_off = SHG::U_EGM2008(tc.r_m, phi, lambda, tc.degree);

        if (!nearly_equal(norm3(g_on), norm3(g_off), kGTol) || !nearly_equal(u_on, u_off, kUTol)) {
            ++failures;
            std::cerr << "[FAIL] workspace toggle mismatch for degree=" << tc.degree
                      << " phi=" << tc.phi_deg << " lambda=" << tc.lambda_deg << std::endl;
            std::cerr << "  |g| on=" << norm3(g_on) << " off=" << norm3(g_off)
                      << " diff=" << std::abs(norm3(g_on) - norm3(g_off)) << std::endl;
            std::cerr << "  U   on=" << u_on << " off=" << u_off
                      << " diff=" << std::abs(u_on - u_off) << std::endl;
        }
    }

    // Performance comparison (informational, non-gating)
    const int iterations = 30;
    const auto perf_case = cases[2]; // medium/high degree representative point
    const double phi = deg2rad(perf_case.phi_deg);
    const double lambda = deg2rad(perf_case.lambda_deg);

    volatile double sink = 0.0;

    SHG::set_workspace_enabled(false);
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        auto g = SHG::g_EGM2008(perf_case.r_m, phi, lambda, perf_case.degree);
        sink += norm3(g);
        sink += SHG::U_EGM2008(perf_case.r_m, phi, lambda, perf_case.degree);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    SHG::set_workspace_enabled(true);
    auto t2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        auto g = SHG::g_EGM2008(perf_case.r_m, phi, lambda, perf_case.degree);
        sink += norm3(g);
        sink += SHG::U_EGM2008(perf_case.r_m, phi, lambda, perf_case.degree);
    }
    auto t3 = std::chrono::high_resolution_clock::now();

    const double off_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double on_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

    std::cout << "Workspace comparison (lower is better): off=" << off_ms
              << " ms, on=" << on_ms << " ms, checksum=" << sink << std::endl;

    SHG::set_workspace_enabled(true); // restore default behavior

    if (failures == 0) {
        std::cout << "Workspace toggle equivalence checks passed." << std::endl;
        return 0;
    }

    std::cerr << failures << " workspace comparison failure(s)." << std::endl;
    return 1;
}
