
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <unistd.h> // for getcwd


#include "SHG.cpp"

int main(int argc, char *argv[])
{
    // print cwd
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
    {
        std::cout << "Current working dir: " << cwd << std::endl;
    }
    else
    {
        std::cerr << "getcwd() error" << std::endl;
    }

    // Example usage
    int l_max = 2;         // Maximum degree
    int m_max = 2;         // Maximum order
    using namespace SHG;

    // read coefficients from binary file
    std::vector<std::vector<double>> C, S;
    if (!read_coefficients_block_binary("EGM2008Coeffs.bin", l_max, m_max, C, S))
    {
        std::cerr << "Failed to read coefficients from binary file." << std::endl;
        return 1;
    }

    double GM = 3986004.415e8 / 1e9; // Earth's gravitational parameter in km^3/s^2
    double a = 6378136.3 / 1e3;      // Earth's reference radius in km
    double f = 298.257223563;        // Earth's flattening factor

    // Compute gravitational potential and acceleration at a sample point
    double r = 7000.0; // Radial distance in km
    double phi=45 * M_PI / 180.0;    // Example geocentric latitude (45 degrees)
    double lambda=0.0 * M_PI / 180.0; // Example longitude
    double V = compute_gravitational_potential(r, phi, lambda, l_max, m_max, C, S, a, GM);
    auto accel = compute_gravitational_acceleration<std::array<double, 3>>(r, phi, lambda, l_max, m_max, C, S, a, GM);

    // Output results
    std::cout.precision(15);
    std::cout << "Gravitational Potential V: " << V << " km^2/s^2" << std::endl;
    std::cout << "Gravitational Acceleration (km/s^2): " << std::endl;
    std::cout << "  a_x: " << accel[0] << std::endl;
    std::cout << "  a_y: " << accel[1] << std::endl;
    std::cout << "  a_z: " << accel[2] << std::endl;    

    return 0;
}