#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "SHG.hpp"

namespace {
constexpr double kPi = 3.14159265358979323846;

double deg2rad(double deg) {
    return deg * (kPi / 180.0);
}

struct Point {
    double r_m;
    double phi_deg;
    double lambda_deg;
};

} // namespace

int main(int argc, char* argv[]) {
    int degree = 1080;
    int iterations = 500;
    std::string coefficient_path;

    if (argc >= 2) {
        degree = std::stoi(argv[1]);
    }
    if (argc >= 3) {
        iterations = std::stoi(argv[2]);
    }
    if (argc >= 4) {
        coefficient_path = argv[3];
    }

    const std::vector<Point> points = {
        {6378137.0, 0.0, 0.0},
        {6378237.0, 38.9072, -77.0369},
        {6378137.0, 45.0, 90.0},
        {6378137.0, -45.0, 0.0},
        {6378137.0, 90.0, 0.0},
        {6388137.0, -33.8688, 151.2093},
        {6378137.0, 35.6762, 139.6503},
        {6478137.0, 10.0, 20.0},
    };

    SHG::set_coefficient_loading_verbose(false);

    // Warm-up load/path validation
    {
        const auto& p = points.front();
        const double phi = deg2rad(p.phi_deg);
        const double lambda = deg2rad(p.lambda_deg);
        (void)SHG::U_EGM2008(p.r_m, phi, lambda, degree, coefficient_path);
    }

    volatile double sink = 0.0;

    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        for (const auto& p : points) {
            const double phi = deg2rad(p.phi_deg);
            const double lambda = deg2rad(p.lambda_deg);
            sink += SHG::U_EGM2008(p.r_m, phi, lambda, degree, coefficient_path);
        }
    }
    const auto t1 = std::chrono::high_resolution_clock::now();

    const double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double calls = static_cast<double>(iterations * static_cast<int>(points.size()));

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "mode=potential, degree=" << degree << ", iterations=" << iterations
              << ", points_per_iter=" << points.size() << std::endl;
    std::cout << "U_EGM2008 total_ms=" << total_ms << ", per_call_ms=" << (total_ms / calls) << std::endl;
    std::cout << "checksum=" << sink << std::endl;

    return 0;
}
