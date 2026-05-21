#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
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

    const std::filesystem::path model_dir = SHG::model_bins_directory();
    const std::filesystem::path model_bin = model_dir / "EGM2008.bin";
    if (!std::filesystem::exists(model_bin)) {
        std::cout << "SKIPPED: EGM2008 artifact bin not available at " << model_bin << std::endl;
        return 0;
    }

    SHG::SHGModel model;
    try {
        model = SHG::get_model("EGM2008");
    } catch (const std::exception& e) {
        std::cerr << "ERROR: failed to load EGM2008 artifact model: " << e.what() << std::endl;
        return 1;
    }

    int failures = 0;
    for (const Case& tc : cases) {
        const double phi = deg2rad(tc.phi_deg);
        const double lambda = deg2rad(tc.lambda_deg);
        const int degree = std::min(tc.degree, model.l_max);
        const int order = std::min(tc.degree, model.m_max);

        const std::array<double, 3> g_model = SHG::acceleration(model, tc.r_m, phi, lambda, degree, order);
        const double u_model = SHG::potential(model, tc.r_m, phi, lambda, degree, order);

        const std::array<double, 3> g_explicit = SHG::g_EGM2008(tc.r_m, phi, lambda, tc.degree);
        const double u_explicit = SHG::U_EGM2008(tc.r_m, phi, lambda, tc.degree);

        const double g_model_mag = norm3(g_model);
        const double g_explicit_mag = norm3(g_explicit);

        const bool accel_ok = nearly_equal_machine(g_model_mag, g_explicit_mag);
        const bool pot_ok = nearly_equal_machine(u_model, u_explicit);

        if (!accel_ok || !pot_ok) {
            ++failures;
            std::cerr << "[FAIL] r=" << tc.r_m << " phi=" << tc.phi_deg << " lambda=" << tc.lambda_deg
                      << " degree=" << tc.degree << std::endl;
            std::cerr << "  accel model=" << g_model_mag << " explicit=" << g_explicit_mag
                      << " diff=" << std::abs(g_model_mag - g_explicit_mag) << std::endl;
            std::cerr << "  pot   model=" << u_model << " explicit=" << u_explicit
                      << " diff=" << std::abs(u_model - u_explicit) << std::endl;
        }
    }

    if (failures == 0) {
        std::cout << "EGM2008 artifact parity checks passed (" << cases.size() << " cases)." << std::endl;
        return 0;
    }

    std::cerr << failures << " parity case(s) failed." << std::endl;
    return 1;
}
