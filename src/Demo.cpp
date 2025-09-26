
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
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

    std::cout << "\n=== Testing EGM2008 Specific Functions ===\n" << std::endl;
    
    using namespace SHG;

    // Test locations
    std::vector<std::array<double, 3>> test_locations = {
        {38.9072, -77.0369, 0.1},    // Washington, D.C.
        {40.7128, -74.0060, 0.0},    // New York City  
        {51.5074, -0.1278, 0.0},     // London
        {35.6762, 139.6503, 0.0},    // Tokyo
        {-33.8688, 151.2093, 0.0}    // Sydney
    };
    
    std::vector<std::string> location_names = {
        "Washington D.C.", "New York City", "London", "Tokyo", "Sydney"
    };
    int deg = 1200;
    std::cout << "=== Testing with EGM2008 model (degree 1200) ===\n" << std::endl;

    std::cout << std::setw(15) << "Location"
              << std::setw(12) << "Lat (°)"
              << std::setw(12) << "Lon (°)"
              << std::setw(12) << "Alt (km)"
              << std::setw(15) << "g_mag (m/s²)"
              << std::setw(18) << "U (m²/s²)"
              << std::endl;
    std::cout << std::string(85, '-') << std::endl;

    for (size_t i = 0; i < test_locations.size(); ++i) {
        double lat = test_locations[i][0];
        double lon = test_locations[i][1]; 
        double alt_km = test_locations[i][2];
        
        // Convert to geocentric coordinates
        double phi = lat * M_PI / 180.0;  // Geocentric latitude in radians
        double lambda = lon * M_PI / 180.0;  // Longitude in radians
        double r = 6378137.0 + alt_km * 1000.0;  // Approximate radius in meters

        try {
            // Use new EGM2008 functions with full degree/order (2190)
            auto g_acc = g_EGM2008(r, phi, lambda, deg);
            double g_mag = std::sqrt(g_acc[0]*g_acc[0] + g_acc[1]*g_acc[1] + g_acc[2]*g_acc[2]);
            double potential = U_EGM2008(r, phi, lambda, deg);
            
            std::cout << std::setw(15) << location_names[i]
                      << std::setw(12) << std::fixed << std::setprecision(4) << lat
                      << std::setw(12) << std::fixed << std::setprecision(4) << lon  
                      << std::setw(12) << std::fixed << std::setprecision(1) << alt_km
                      << std::setw(15) << std::fixed << std::setprecision(6) << g_mag
                      << std::setw(18) << std::fixed << std::setprecision(1) << potential
                      << std::endl;
                      
        } catch (const std::exception& e) {
            std::cout << std::setw(15) << location_names[i] << " ERROR: " << e.what() << std::endl;
        }
    }
    
    std::cout << "\n=== Testing with different degrees at equator (100m altitude) ===\n" << std::endl;
    
    // Test with different degrees at same location
    double r = 6378137.0 + 100.0;  // 100m altitude
    double phi = 0.0;  // Equator
    double lambda = 0.0;  // Prime meridian
    
    std::vector<int> test_degrees = {2, 10, 50, 180, 360, 1800};
    
    std::cout << std::setw(10) << "Degree"
              << std::setw(18) << "g_magnitude (m/s²)"
              << std::setw(20) << "U (×10¹² m²/s²)"
              << std::endl;
    std::cout << std::string(75, '-') << std::endl;
    
    for (int degree : test_degrees) {
        try {
            auto g_acc = g_EGM2008(r, phi, lambda, degree);
            double g_mag = std::sqrt(g_acc[0]*g_acc[0] + g_acc[1]*g_acc[1] + g_acc[2]*g_acc[2]);
            double potential = U_EGM2008(r, phi, lambda, degree);
            
            std::cout << std::setw(10) << degree
                      << std::setw(18) << std::fixed << std::setprecision(8) << g_mag
                      << std::setw(20) << std::fixed << std::setprecision(6) << potential/1e12
                      << std::setw(15);
            
            std::cout << std::endl;
            
        } catch (const std::exception& e) {
            std::cout << std::setw(10) << degree << " ERROR: " << e.what() << std::endl;
        }
    }
    
    std::cout << "\n=== Testing with custom coefficient path ===\n" << std::endl;
    std::cout << "Testing explicit binary file path:" << std::endl;
    
    try {
        auto g_acc = g_EGM2008(r, phi, lambda, 100, "EGM2008Coeffs.bin");
        double g_mag = std::sqrt(g_acc[0]*g_acc[0] + g_acc[1]*g_acc[1] + g_acc[2]*g_acc[2]);
        double potential = U_EGM2008(r, phi, lambda, 100, "EGM2008Coeffs.bin");
        
        std::cout << "With degree 100:" << std::endl;
        std::cout << "g = [" << g_acc[0] << ", " << g_acc[1] << ", " << g_acc[2] << "] m/s²" << std::endl;
        std::cout << "|g| = " << g_mag << " m/s²" << std::endl;  
        std::cout << "U = " << potential << " m²/s²" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
    }

    return 0;
}