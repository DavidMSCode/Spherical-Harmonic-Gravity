#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "SHG.hpp"

namespace {
struct RegressionCase {
    std::string name;
    double r_m;
    double phi_deg;
    double lambda_deg;
    int degree;
    double expected_accel_mps2;
    double expected_potential_m2ps2;
};

constexpr double kPi = 3.14159265358979323846;
constexpr double kAccelTol = 1e-12;
constexpr double kPotTol = 1e-6;

double deg2rad(double deg) {
    return deg * (kPi / 180.0);
}

bool nearly_equal(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

bool parse_case(const std::string& line, RegressionCase& out) {
    std::stringstream ss(line);
    std::string token;
    std::vector<std::string> cols;

    while (std::getline(ss, token, ',')) {
        cols.push_back(token);
    }

    if (cols.size() != 7) {
        return false;
    }

    out.name = cols[0];
    out.r_m = std::stod(cols[1]);
    out.phi_deg = std::stod(cols[2]);
    out.lambda_deg = std::stod(cols[3]);
    out.degree = std::stoi(cols[4]);
    out.expected_accel_mps2 = std::stod(cols[5]);
    out.expected_potential_m2ps2 = std::stod(cols[6]);
    return true;
}

std::vector<RegressionCase> load_cases_from_csv(const std::string& csv_path) {
    std::ifstream in(csv_path);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open baseline CSV: " + csv_path);
    }

    std::vector<RegressionCase> cases;
    std::string line;
    bool first_line = true;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (first_line) {
            first_line = false;
            continue; // header
        }

        RegressionCase rc{};
        if (!parse_case(line, rc)) {
            throw std::runtime_error("Invalid CSV row: " + line);
        }
        cases.push_back(rc);
    }

    if (cases.empty()) {
        throw std::runtime_error("Baseline CSV contains no regression rows: " + csv_path);
    }

    return cases;
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <baseline_csv_path> [coefficient_path]" << std::endl;
        return 2;
    }

    const std::string baseline_csv_path = argv[1];
    const std::string coefficient_path = (argc > 2) ? argv[2] : "";

    std::vector<RegressionCase> cases;
    try {
        cases = load_cases_from_csv(baseline_csv_path);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 2;
    }

    int failures = 0;

    for (const auto& tc : cases) {
        const double phi = deg2rad(tc.phi_deg);
        const double lambda = deg2rad(tc.lambda_deg);

        const auto g = SHG::g_EGM2008(tc.r_m, phi, lambda, tc.degree, coefficient_path);
        const double g_mag = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
        const double potential = SHG::U_EGM2008(tc.r_m, phi, lambda, tc.degree, coefficient_path);

        const bool accel_ok = nearly_equal(g_mag, tc.expected_accel_mps2, kAccelTol);
        const bool pot_ok = nearly_equal(potential, tc.expected_potential_m2ps2, kPotTol);

        if (!accel_ok || !pot_ok) {
            ++failures;
            std::cerr << "[FAIL] " << tc.name << '\n';
            if (!accel_ok) {
                std::cerr << "  acceleration expected=" << tc.expected_accel_mps2
                          << " got=" << g_mag
                          << " |diff|=" << std::abs(g_mag - tc.expected_accel_mps2)
                          << " tol=" << kAccelTol << '\n';
            }
            if (!pot_ok) {
                std::cerr << "  potential    expected=" << tc.expected_potential_m2ps2
                          << " got=" << potential
                          << " |diff|=" << std::abs(potential - tc.expected_potential_m2ps2)
                          << " tol=" << kPotTol << '\n';
            }
        }
    }

    if (failures == 0) {
        std::cout << "All EGM2008 regression checks passed (" << cases.size() << " cases)." << std::endl;
        return 0;
    }

    std::cerr << failures << " regression case(s) failed." << std::endl;
    return 1;
}
