#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "SHG.hpp"

namespace {
struct BaselineInput {
    std::string name;
    double r_m;
    double phi_deg;
    double lambda_deg;
    int degree;
};

constexpr double kPi = 3.14159265358979323846;

double deg2rad(double deg) {
    return deg * (kPi / 180.0);
}

} // namespace

int main(int argc, char* argv[]) {
    const std::string output_csv_path = (argc > 1) ? argv[1] : "";
    const std::string coefficient_path = (argc > 2) ? argv[2] : "";

    const std::vector<BaselineInput> cases = {
        {"equator_sea_level_d360", 6378137.0, 0.0, 0.0, 360},
        {"washington_dc_d360", 6378237.0, 38.9072, -77.0369, 360},
        {"lat45_lon90_d180", 6378137.0, 45.0, 90.0, 180},
        {"south45_lon0_d180", 6378137.0, -45.0, 0.0, 180},
        {"north_pole_d180", 6378137.0, 90.0, 0.0, 180},
        {"equator_plus100m_d720", 6378237.0, 0.0, 0.0, 720},
        {"tokyo_sea_level_d1080", 6378137.0, 35.6762, 139.6503, 1080},
        {"sydney_10km_d720", 6388137.0, -33.8688, 151.2093, 720},
    };

    std::ofstream file_out;
    std::ostream* out = &std::cout;

    if (!output_csv_path.empty()) {
        file_out.open(output_csv_path);
        if (!file_out.is_open()) {
            std::cerr << "ERROR: Could not open output file: " << output_csv_path << std::endl;
            return 1;
        }
        out = &file_out;
    }

    *out << "name,r_m,phi_deg,lambda_deg,degree,acceleration_mps2,potential_m2ps2\n";
    out->setf(std::ios::fmtflags(0), std::ios::floatfield);
    *out << std::setprecision(std::numeric_limits<double>::max_digits10);

    for (const auto& tc : cases) {
        const double phi = deg2rad(tc.phi_deg);
        const double lambda = deg2rad(tc.lambda_deg);

        const auto g = SHG::g_EGM2008(tc.r_m, phi, lambda, tc.degree, coefficient_path);
        const double g_mag = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
        const double potential = SHG::U_EGM2008(tc.r_m, phi, lambda, tc.degree, coefficient_path);

        *out << tc.name << ','
             << tc.r_m << ','
             << tc.phi_deg << ','
             << tc.lambda_deg << ','
             << tc.degree << ','
             << g_mag << ','
             << potential << '\n';
    }

    if (!output_csv_path.empty()) {
        std::cout << "Wrote full-precision baseline CSV to " << output_csv_path << std::endl;
    }

    return 0;
}
