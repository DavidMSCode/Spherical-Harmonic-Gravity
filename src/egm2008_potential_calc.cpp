/*
 * EGM2008 Gravitational Potential Calculator
 * -----------------------------------------
 * Interactive program to compute gravitational potential using EGM2008 model.
 * 
 * Author: DavidMSCode
 * Project: Spherical-Harmonic-Gravity
 * License: MIT
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "SHG.hpp"

void printHeader() {
    std::cout << "===============================================" << std::endl;
    std::cout << "  EGM2008 Gravitational Potential Calculator" << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << "Computes gravitational potential using the" << std::endl;
    std::cout << "EGM2008 Earth Gravitational Model." << std::endl;
    std::cout << std::endl;
}

void printUsage() {
    std::cout << "Input formats:" << std::endl;
    std::cout << "  Radius:     geocentric radius in meters (e.g., 6378237 for 100m altitude)" << std::endl;
    std::cout << "  Latitude:   geocentric latitude in degrees (e.g., 38.9072)" << std::endl;
    std::cout << "  Longitude:  longitude in degrees (e.g., -77.0369)" << std::endl;
    std::cout << "  Max Degree: 2-2190 (recommend ≤1900 for stability)" << std::endl;
    std::cout << std::endl;
}

bool getInput(double& r, double& phi, double& lambda, int& max_degree, std::string& coeff_path) {
    std::cout << "Enter geocentric radius (meters): ";
    if (!(std::cin >> r)) {
        std::cin.clear();
        std::cin.ignore(10000, '\n');
        return false;
    }
    
    std::cout << "Enter geocentric latitude (degrees): ";
    double lat_deg;
    if (!(std::cin >> lat_deg)) {
        std::cin.clear();
        std::cin.ignore(10000, '\n');
        return false;
    }
    phi = lat_deg * M_PI / 180.0; // Convert to radians
    
    std::cout << "Enter longitude (degrees): ";
    double lon_deg;
    if (!(std::cin >> lon_deg)) {
        std::cin.clear();
        std::cin.ignore(10000, '\n');
        return false;
    }
    lambda = lon_deg * M_PI / 180.0; // Convert to radians
    
    std::cout << "Enter max degree (2-2190, recommend ≤1900): ";
    if (!(std::cin >> max_degree)) {
        std::cin.clear();
        std::cin.ignore(10000, '\n');
        return false;
    }
    
    std::cout << "Enter coefficient file path (or press Enter for default): ";
    std::cin.ignore(); // Clear the newline from previous input
    std::getline(std::cin, coeff_path);
    
    // Validate inputs
    if (lat_deg < -90 || lat_deg > 90) {
        std::cout << "ERROR: Latitude must be between -90 and 90 degrees." << std::endl;
        return false;
    }
    
    if (lon_deg < -180 || lon_deg > 180) {
        std::cout << "ERROR: Longitude must be between -180 and 180 degrees." << std::endl;
        return false;
    }
    
    if (max_degree < 2 || max_degree > 2190) {
        std::cout << "ERROR: Max degree must be between 2 and 2190." << std::endl;
        return false;
    }
    
    return true;
}

void printResults(double potential, double r, double phi, double lambda, int max_degree, bool verbose = true) {
    if (verbose) {
        std::cout << std::endl;
        std::cout << "========== RESULTS ==========" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Geocentric radius: " << r << " m" << std::endl;
        std::cout << "Geocentric latitude: " << phi * 180.0 / M_PI << "°" << std::endl;
        std::cout << "Longitude: " << lambda * 180.0 / M_PI << "°" << std::endl;
        std::cout << "EGM2008 degree/order: " << max_degree << std::endl;
        std::cout << std::endl;
        std::cout << "Gravitational potential: " << std::setw(15) << std::setprecision(1) 
                  << potential << " m²/s²" << std::endl;
        
        if (max_degree > 1900) {
            std::cout << std::endl;
            std::cout << "NOTE: Using degree > 1900 may have numerical instability." << std::endl;
        }
    } else {
        // Command-line output: just the potential value, nicely formatted
        std::cout << std::fixed << std::setprecision(1) << potential << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // Check for help flag
    if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        std::cout << "EGM2008 Gravitational Potential Calculator" << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "  Interactive mode:     " << argv[0] << std::endl;
        std::cout << "  Command-line mode:    " << argv[0] << " [--verbose|-v] <r> <phi> <lambda> <degree> [coeff_path]" << std::endl;
        std::cout << std::endl;
        std::cout << "Parameters:" << std::endl;
        std::cout << "  --verbose, -v Show coefficient loading messages" << std::endl;
        std::cout << "  r             Geocentric radius in meters" << std::endl;
        std::cout << "  phi           Geocentric latitude in degrees" << std::endl;
        std::cout << "  lambda        Longitude in degrees" << std::endl;
        std::cout << "  degree        Maximum degree/order to use" << std::endl;
        std::cout << "  coeff_path    Optional path to coefficient file" << std::endl;
        std::cout << std::endl;
        std::cout << "Output (command-line mode):" << std::endl;
        std::cout << "  Single value: gravitational potential in m^2/s^2" << std::endl;
        return 0;
    }
    
    // Check for verbose flag and command-line mode
    bool verbose_flag = false;
    int arg_offset = 1;
    
    if (argc >= 2 && (std::string(argv[1]) == "--verbose" || std::string(argv[1]) == "-v")) {
        verbose_flag = true;
        arg_offset = 2;
    }
    
    if (argc >= arg_offset + 4) {
        double r = std::atof(argv[arg_offset]);
        double phi_deg = std::atof(argv[arg_offset + 1]);
        double lambda_deg = std::atof(argv[arg_offset + 2]);
        int max_degree = std::atoi(argv[arg_offset + 3]);
        std::string coeff_path = (argc >= arg_offset + 5) ? argv[arg_offset + 4] : "";
        
        // Convert to radians
        double phi = phi_deg * M_PI / 180.0;
        double lambda = lambda_deg * M_PI / 180.0;
        
        // Set verbosity based on flag
        SHG::set_coefficient_loading_verbose(verbose_flag);
        
        try {
            double potential = SHG::U_EGM2008(r, phi, lambda, max_degree, coeff_path);
            
            // Compact output for command-line mode
            printResults(potential, r, phi, lambda, max_degree, false);
            
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
        
        return 0;
    }
    
    // Interactive mode
    printHeader();
    printUsage();
    
    // Enable verbose coefficient loading for interactive mode
    SHG::set_coefficient_loading_verbose(true);
    
    double r, phi, lambda;
    int max_degree;
    std::string coeff_path;
    char continue_calc = 'y';
    
    while (continue_calc == 'y' || continue_calc == 'Y') {
        try {
            if (!getInput(r, phi, lambda, max_degree, coeff_path)) {
                std::cout << "Invalid input. Please try again." << std::endl;
                continue;
            }
            
            std::cout << std::endl << "Computing gravitational potential..." << std::endl;
            
            // Compute potential
            double potential = SHG::U_EGM2008(r, phi, lambda, max_degree, coeff_path);
            
            // Verbose output for interactive mode
            printResults(potential, r, phi, lambda, max_degree, true);
            
        } catch (const std::exception& e) {
            std::cout << "ERROR: " << e.what() << std::endl;
        }
        
        std::cout << std::endl;
        std::cout << "Calculate another location? (y/n): ";
        std::cin >> continue_calc;
    }
    
    std::cout << std::endl << "Thank you for using the EGM2008 Potential Calculator!" << std::endl;
    return 0;
}