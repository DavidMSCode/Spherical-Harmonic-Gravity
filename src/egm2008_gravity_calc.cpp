/*
 * EGM2008 Gravitational Acceleration Calculator
 * --------------------------------------------
 * Interactive program to compute gravitational acceleration using EGM2008 model.
 * 
 * Author: DavidMSCode
 * Project: Spherical-Harmonic-Gravity
 * License: MIT
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <array>
#include "SHG.hpp"

void printHeader() {
    std::cout << "===============================================" << std::endl;
    std::cout << "  EGM2008 Gravitational Acceleration Calculator" << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << "Computes gravitational acceleration using the" << std::endl;
    std::cout << "EGM2008 Earth Gravitational Model." << std::endl;
    std::cout << std::endl;
}

void printUsage() {
    std::cout << "Input formats:" << std::endl;
    std::cout << "  Radius:         geocentric radius in meters (e.g., 6378237 for sea level)" << std::endl;
    std::cout << "  Geo. Latitude:  geocentric latitude in degrees (e.g., 38.9072 for Washington DC)" << std::endl;
    std::cout << "  Longitude:      longitude in degrees (e.g., -77.0369 for Washington DC)" << std::endl;
    std::cout << "  Max Degree:     2-2190 (recommend ≤1900 for stability)" << std::endl;
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
    if (!(std::cin >> phi)) {
        std::cin.clear();
        std::cin.ignore(10000, '\n');
        return false;
    }
    
    std::cout << "Enter longitude (degrees): ";
    if (!(std::cin >> lambda)) {
        std::cin.clear();
        std::cin.ignore(10000, '\n');
        return false;
    }
    
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
    if (r <= 0) {
        std::cout << "ERROR: Radius must be positive." << std::endl;
        return false;
    }
    
    if (phi < -90 || phi > 90) {
        std::cout << "ERROR: Geocentric latitude must be between -90 and 90 degrees." << std::endl;
        return false;
    }
    
    if (lambda < -180 || lambda > 180) {
        std::cout << "ERROR: Longitude must be between -180 and 180 degrees." << std::endl;
        return false;
    }
    
    if (max_degree < 2 || max_degree > 2190) {
        std::cout << "ERROR: Max degree must be between 2 and 2190." << std::endl;
        return false;
    }
    
    return true;
}

void printResults(const std::array<double, 3>& g_acc, double r, double phi, double lambda, int max_degree, bool verbose = true) {
    double g_mag = std::sqrt(g_acc[0]*g_acc[0] + g_acc[1]*g_acc[1] + g_acc[2]*g_acc[2]);
    
    if (verbose) {
        std::cout << std::endl;
        std::cout << "========== RESULTS ==========" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Location: (r=" << r << "m, φ=" << phi << "°, λ=" << lambda << "°)" << std::endl;
        std::cout << "EGM2008 degree/order: " << max_degree << std::endl;
        std::cout << std::endl;
        std::cout << "Gravitational acceleration components:" << std::endl;
        std::cout << "  gI (radial):      " << std::setw(12) << g_acc[0] << " m/s²" << std::endl;
        std::cout << "  gJ (south):       " << std::setw(12) << g_acc[1] << " m/s²" << std::endl;
        std::cout << "  gK (east):        " << std::setw(12) << g_acc[2] << " m/s²" << std::endl;
        std::cout << std::endl;
        std::cout << "Magnitude:          " << std::setw(12) << g_mag << " m/s²" << std::endl;
        
        if (max_degree > 1900) {
            std::cout << "NOTE: Using degree > 1900 may have numerical instability." << std::endl;
        }
        
        double theoretical_g = 9.80665; // Standard gravity
        double diff_percent = ((g_mag - theoretical_g) / theoretical_g) * 100.0;
        std::cout << "Difference from standard gravity (9.80665 m/s²): " 
                  << std::setprecision(3) << diff_percent << "%" << std::endl;
    } else {
        // Command-line output: just the gravity magnitude, nicely formatted
        std::cout << std::fixed << std::setprecision(6) << g_mag << std::endl;
    }
}

void printCmdUsage() {
    std::cout << "EGM2008 Gravitational Acceleration Calculator" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage:" << std::endl;
    std::cout << "  Interactive mode:     ./GravityCalc" << std::endl;
    std::cout << "  Command-line mode:    ./GravityCalc [--verbose|-v] <r> <phi> <lambda> <degree> [coeff_path]" << std::endl;
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
    std::cout << "  Single value: gravitational acceleration magnitude in m/s^2" << std::endl;
}

int main(int argc, char* argv[]) {
    // Check for help flag
    if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        printCmdUsage();
        return 0;
    }
    
    // Check for command line arguments (with or without verbose flag)
    bool verbose_flag = false;
    int arg_offset = 1;
    
    if (argc >= 2 && (std::string(argv[1]) == "--verbose" || std::string(argv[1]) == "-v")) {
        verbose_flag = true;
        arg_offset = 2;
    }
    
    if (argc >= arg_offset + 4) {
        // Command line mode
        try {
            double r = std::stod(argv[arg_offset]);
            double phi_deg = std::stod(argv[arg_offset + 1]);
            double lambda_deg = std::stod(argv[arg_offset + 2]);
            int max_degree = std::stoi(argv[arg_offset + 3]);
            std::string coeff_path = (argc >= arg_offset + 5) ? argv[arg_offset + 4] : "";
            
            // Validate inputs
            if (r <= 0) {
                std::cerr << "ERROR: Radius must be positive." << std::endl;
                return 1;
            }
            if (phi_deg < -90 || phi_deg > 90) {
                std::cerr << "ERROR: Geocentric latitude must be between -90 and 90 degrees." << std::endl;
                return 1;
            }
            if (lambda_deg < -180 || lambda_deg > 180) {
                std::cerr << "ERROR: Longitude must be between -180 and 180 degrees." << std::endl;
                return 1;
            }
            if (max_degree < 2 || max_degree > 2190) {
                std::cerr << "ERROR: Max degree must be between 2 and 2190." << std::endl;
                return 1;
            }
            
            // Convert to radians for computation
            double phi = phi_deg * M_PI / 180.0;
            double lambda = lambda_deg * M_PI / 180.0;
            
            // Set verbosity based on flag
            SHG::set_coefficient_loading_verbose(verbose_flag);
            
            // Compute acceleration
            auto g_acc = SHG::g_EGM2008(r, phi, lambda, max_degree, coeff_path);
            
            // Print compact results
            printResults(g_acc, r, phi_deg, lambda_deg, max_degree, false);
            
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return 1;
        }
        
        return 0;
    }
    
    // Interactive mode
    if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        printCmdUsage();
        return 0;
    }
    
    printHeader();
    printUsage();
    
    // Enable verbose coefficient loading for interactive mode  
    SHG::set_coefficient_loading_verbose(true);
    
    double r, phi_deg, lambda_deg;
    int max_degree;
    std::string coeff_path;
    char continue_calc = 'y';
    
    while (continue_calc == 'y' || continue_calc == 'Y') {
        try {
            if (!getInput(r, phi_deg, lambda_deg, max_degree, coeff_path)) {
                std::cout << "Invalid input. Please try again." << std::endl;
                continue;
            }
            
            // Convert to radians for computation
            double phi = phi_deg * M_PI / 180.0;
            double lambda = lambda_deg * M_PI / 180.0;
            
            std::cout << std::endl << "Computing gravitational acceleration..." << std::endl;
            
            // Compute acceleration
            auto g_acc = SHG::g_EGM2008(r, phi, lambda, max_degree, coeff_path);
            
            printResults(g_acc, r, phi_deg, lambda_deg, max_degree, true);
            
        } catch (const std::exception& e) {
            std::cout << "ERROR: " << e.what() << std::endl;
        }
        
        std::cout << std::endl;
        std::cout << "Calculate another location? (y/n): ";
        std::cin >> continue_calc;
    }
    
    std::cout << std::endl << "Thank you for using the EGM2008 Gravity Calculator!" << std::endl;
    return 0;
}