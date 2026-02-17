/*
 * Spherical Harmonic Gravity (SHG.hpp)
 * -------------------------------------
 * Computes gravitational potential and acceleration using normalized spherical harmonic coefficients.
 *
 * Project: Spherical-Harmonic-Gravity
 * License: MIT (see LICENSE file for details)
 * Author: DavidMSCode
 * Created: 2025
 *
 * This file is part of a permissively licensed project. You may use, modify, and distribute
 * this code for personal or commercial purposes, provided the license is included with any
 * modified source code. No warranty is provided; see LICENSE for details.
 *
 */


#ifndef SHG_HPP
#define SHG_HPP

#include <vector>
#include <string>
#include <utility> // for std::pair

namespace SHG {
    
    // Physical constants and reference parameters
    
    // WGS84 Constants
    constexpr double WGS84_A = 6378137.0;           // WGS84 semi-major axis in meters
    constexpr double WGS84_A_KM = 6378.137;         // WGS84 semi-major axis in kilometers
    constexpr double WGS84_GM = 3.986004418e14;     // WGS84 gravitational parameter in m³/s²
    constexpr double WGS84_GM_KM = 3.986004418e5;   // WGS84 gravitational parameter in km³/s²
    constexpr double WGS84_F = 1.0 / 298.257223563; // WGS84 flattening
    
    // EGM2008 Constants (from README_WGS84_2.pdf)
    constexpr double EGM2008_A = 6378136.3;         // EGM2008 semi-major axis in meters
    constexpr double EGM2008_A_KM = 6378.1363;      // EGM2008 semi-major axis in kilometers
    constexpr double EGM2008_GM = 3.986004415e14;   // EGM2008 gravitational parameter in m³/s²
    constexpr double EGM2008_GM_KM = 3.986004415e5; // EGM2008 gravitational parameter in km³/s²
    constexpr int EGM2008_MAX_DEGREE = 2190;        // EGM2008 maximum degree
    constexpr int EGM2008_MAX_ORDER = 2190;         // EGM2008 maximum order
    
    /**
     * @brief Reads a block (starting from 0,0) of spherical harmonic coefficients from a binary file.
     * @param filename The path to the binary file containing the coefficients.
     * @param l Maximum degree of the coefficients to read.
     * @param m Maximum order of the coefficients to read.
     * @param C_block Output 2D vector to store cosine coefficients C[l][m].
     * @param S_block Output 2D vector to store sine coefficients S[l][m].
     * @return True if the coefficients were read successfully, false otherwise.
     */
	bool read_coefficients_block_binary(const std::string &filename, int l, int m, std::vector<std::vector<double>> &C_block, std::vector<std::vector<double>> &S_block);

    /**
     * @brief Reads a specific spherical harmonic coefficient (C and S) from a binary file.
     * @param filename The path to the binary file containing the coefficients.
     * @param l Degree of the coefficient to read.
     * @param m Order of the coefficient to read.
     * @param C_val Output variable to store the cosine coefficient C[l][m].
     * @param S_val Output variable to store the sine coefficient S[l][m].
     * @return True if the coefficient was read successfully, false otherwise.
     */
    bool read_coefficient_binary(const std::string &filename, int l, int m, double &C_val, double &S_val);

    /**
     * @brief Writes spherical harmonic coefficients to a binary file for efficient storage.
     * @param filename The path to the binary file to write the coefficients to.
     * @param C 2D vector containing cosine coefficients C[l][m].
     * @param S 2D vector containing sine coefficients S[l][m].
     * @param l_max Maximum degree of the coefficients.
     * @param m_max Maximum order of the coefficients.
     * @return True if the coefficients were written successfully, false otherwise.
     */
    bool write_coefficients_binary(const std::string &filename, const std::vector<std::vector<double>> &C, const std::vector<std::vector<double>> &S, int l_max, int m_max);

    /**
     * @brief Computes normalized associated Legendre polynomials P(l,m,sin(phi)) for all degrees and orders up to l_max and m_max at geocentric latitude phi.
     * @param l_max Maximum degree.
     * @param m_max Maximum order.
     * @param phi Geocentric latitude in radians.
     * @return 2D vector P where P[l][m] corresponds to degree l and order m.
     * 
     * Note: This is equivalent to P(l,m,cos(theta)) where theta is the colatitude.
     * Because this uses the unmodified forward column method, the polynomials may become numerically unstable for m > ~1900.
     */
    std::vector<std::vector<double>> Plm_bar(int l_max, int m_max, double phi);

    /**
     * @brief Computes vector of m*Tan(phi) for m=0 to m_max using recursion formulation.
     * @param m_max Maximum order.
     * @param phi Geocentric latitude in radians.
     * @return Vector of size m_max+1 where the index corresponds to m.
     */
    std::vector<double> recursive_tangent(int m_max, double phi);

    /**
     * @brief Computes vectors of sin(m*lambda) and cos(m*lambda) for m=0 to m_max using recursion.
     * @param m_max Maximum order.
     * @param lambda Longitude in radians.
     * @return Pair of vectors: first is sin(m*lambda), second is cos(m*lambda).
     */
    std::pair<std::vector<double>, std::vector<double>> recursive_sine_cosine(int m_max, double lambda);

    /**
     * @brief Reads EGM2008 spherical harmonic coefficients from a text file.
     * @param filename The path to the text file containing the coefficients.
     * @param l_max Maximum degree of the coefficients to read.
     * @param m_max Maximum order of the coefficients to read.
     * @param C Output 2D vector to store cosine coefficients C[l][m].
     * @param S Output 2D vector to store sine coefficients S[l][m].
     * @return True if the coefficients were read successfully, false otherwise.
     *
     * The text file should have the following format:
     * l m C(l,m) S(l,m) {RMS error C} {RMS error S}
     * 2    0   -0.484165143790815D-03   0.000000000000000D+00    0.7481239490D-11    0.0000000000D+00
     *
     * where l is the degree, m is the order, C(l,m) is the cosine coefficient,
     * S(l,m) is the sine coefficient, and RMS error C and RMS error S are
     * optional and can be ignored. The function returns true if the file was
     * loaded successfully, and false otherwise. The file is stored ordered by
     * degree then order, so all coefficients for degree 0 are listed first,
     * followed by all coefficients for degree 1, and so on. Within each degree,
     * the coefficients are listed in order of increasing order.
     */
    bool read_EGM2008_coefficients_text(const std::string &filename, int l_max, int m_max, std::vector<std::vector<double>> &C, std::vector<std::vector<double>> &S);
    
    template <typename T>
    T g(double r, double phi, double lambda, int l_max, int m_max, const std::vector<std::vector<double>> &C, const std::vector<std::vector<double>> &S, double a, double GM);

    double pot(double r, double phi, double lambda, int l_max, int m_max, const std::vector<std::vector<double>> &C, const std::vector<std::vector<double>> &S, double a, double GM);

    std::array<double, 3> cartesian_to_geocentric(double x, double y, double z);

    template <typename T>
    T gravitational_acceleration_from_cartesian(double x, double y, double z, int l_max, int m_max, const std::vector<std::vector<double>> &C, const std::vector<std::vector<double>> &S, double a, double GM);

    // EGM2008-specific functions
    /**
     * @brief Computes gravitational acceleration using EGM2008 model with user-specified degree.
     * Automatically loads coefficients from binary or text file with fallback.
     * @param r Geocentric radius in meters.
     * @param phi Geocentric latitude in radians.
     * @param lambda Longitude in radians.
     * @param max_degree Maximum degree (and order) to use. Must be <= 2190.
     * @param coefficient_path Optional path to EGM2008 coefficient file (.bin or .txt). 
     *                        If empty, searches for "EGM2008Coeffs.bin" then "EGM2008_to2190_TideFree.txt"
     * @return Array of gravitational acceleration components [aI, aJ, aK] in m/s².
     */
    std::array<double, 3> g_EGM2008(double r, double phi, double lambda, int max_degree, const std::string& coefficient_path = "");

    /**
     * @brief Computes gravitational acceleration using EGM2008 model with user-specified degree.
     * Automatically loads coefficients from binary or text file with fallback.
     * @param r Geocentric radius in kilometers.
     * @param phi Geocentric latitude in radians.
     * @param lambda Longitude in radians.
     * @param max_degree Maximum degree (and order) to use. Must be <= 2190.
     * @param coefficient_path Optional path to EGM2008 coefficient file (.bin or .txt). 
     *                        If empty, searches for "EGM2008Coeffs.bin" then "EGM2008_to2190_TideFree.txt"
     * @return Array of gravitational acceleration components [aI, aJ, aK] in km/s².
     */
    std::array<double, 3> g_EGM2008_KM(double r, double phi, double lambda, int max_degree, const std::string &coefficient_path="");

    /**
     * @brief Computes gravitational potential using EGM2008 model with user-specified degree.
     * Automatically loads coefficients from binary or text file with fallback.
     * @param r Geocentric radius in meters.
     * @param phi Geocentric latitude in radians.
     * @param lambda Longitude in radians.
     * @param max_degree Maximum degree (and order) to use. Must be <= 2190.
     * @param coefficient_path Optional path to EGM2008 coefficient file (.bin or .txt).
     *                        If empty, searches for "EGM2008Coeffs.bin" then "EGM2008_to2190_TideFree.txt"
     * @return Gravitational potential in m²/s².
     */
    double U_EGM2008(double r, double phi, double lambda, int max_degree, const std::string& coefficient_path = "");

    /**
     * @brief Computes gravitational potential using EGM2008 model with user-specified degree.
     * Automatically loads coefficients from binary or text file with fallback.
     * @param r Geocentric radius in kilometers.
     * @param phi Geocentric latitude in radians.
     * @param lambda Longitude in radians.
     * @param max_degree Maximum degree (and order) to use. Must be <= 2190.
     * @param coefficient_path Optional path to EGM2008 coefficient file (.bin or .txt).
     *                        If empty, searches for "EGM2008Coeffs.bin" then "EGM2008_to2190_TideFree.txt"
     * @return Gravitational potential in km²/s².
     */
    double U_EGM2008_KM(double r, double phi, double lambda, int max_degree, const std::string &coefficient_path = "");

    /**
     * @brief Loads EGM2008 coefficients with automatic fallback from binary to text format.
     * If binary file doesn't exist, attempts to load from text file and create binary.
     * @param coefficient_path Path to coefficient file (.bin or .txt). If empty, uses default names.
     * @param C Output 2D vector to store cosine coefficients C[l][m].
     * @param S Output 2D vector to store sine coefficients S[l][m].
     * @return True if coefficients were loaded successfully, false otherwise.
     */
    bool load_EGM2008_coefficients(const std::string& coefficient_path, std::vector<std::vector<double>>& C, std::vector<std::vector<double>>& S);

    /**
     * @brief Controls verbosity of coefficient loading messages.
     * @param verbose If true (default), shows coefficient loading messages. If false, suppresses them.
     */
    void set_coefficient_loading_verbose(bool verbose);
}

#include "SHG.tpp" // Include template implementations

#endif // SHG_HPP