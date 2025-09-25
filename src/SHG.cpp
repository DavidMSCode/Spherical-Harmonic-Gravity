/*
 * Spherical Harmonic Gravity (SHG.cpp)
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
 * Code derived from the following sources:
 * Spherical Harmonic Gravity implementation by David A. Vallado
 * - "Fundamentals of Astrodynamics and Applications" by David A. Vallado, 4th Edition
 *
 * Standard forward column method for computing Normalized Associated Legendre Functions (ALFs)
 *  - S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw summation and
 *    the recursive computation of very high degree and order normalized associated Legendre
 *    functions", Journal of Geodesy, vol. 76, pp. 279–299, 2002. https://doi.org/10.1007/s00190-002-0245-x
 *
 * Note that is not using the more stable modified methods described in Holmes and Featherstone
 *  for simplicity and the ALFs are not numerically stable for M>1900 in this method.
 *
 * Don't be confused by the colatitude (theta) vs latitude (phi) notation in Vallado's book. The
 * ALFs of the sine of the geocentric latitude (phi) are the exact same as the ALFs of the cosine
 * of the colatitude (theta).
 *
 * Nomenclature follows Vallado's conventions where applicable.
 * l - degree
 * m - order
 * C - cosine coefficient
 * S - sine coefficient
 * P - associated Legendre polynomial
 * r - radial distance
 * phi - geocentric latitude
 * lambda - longitude
 * GM - gravitational parameter
 * a - reference radius
 */

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility> // for std::pair
#include <chrono>
#include <unistd.h> // for getcwd

#include "SHG.hpp"

namespace SHG
{
    // Reads a block of coefficients from (0,0) to (l,m) from the binary file
    // Returns true if successful, false otherwise
    bool read_coefficients_block_binary(const std::string &filename, int l, int m, std::vector<std::vector<double>> &C_block, std::vector<std::vector<double>> &S_block)
    {
        std::ifstream in(filename, std::ios::binary);
        if (!in.is_open())
            return false;
        int l_max = 0, m_max = 0;
        in.read(reinterpret_cast<char *>(&l_max), sizeof(int));
        in.read(reinterpret_cast<char *>(&m_max), sizeof(int));
        if (l > l_max || m > m_max)
        {
            in.close();
            return false;
        }
        C_block.resize(l + 1, std::vector<double>(m + 1, 0.0));
        S_block.resize(l + 1, std::vector<double>(m + 1, 0.0));
        // Efficiently read C coefficients block
        size_t full_row = m_max + 1;
        size_t block_rows = l + 1;
        size_t block_cols = m + 1;
        std::vector<double> C_buffer(block_rows * full_row);
        std::streamoff offset_C_start = sizeof(int) * 2;
        in.seekg(offset_C_start, std::ios::beg);
        in.read(reinterpret_cast<char *>(C_buffer.data()), sizeof(double) * block_rows * full_row);
        for (int i = 0; i < block_rows; ++i)
        {
            for (int j = 0; j < block_cols; ++j)
            {
                C_block[i][j] = C_buffer[i * full_row + j];
            }
        }
        // Efficiently read S coefficients block
        std::vector<double> S_buffer(block_rows * full_row);
        std::streamoff offset_S_start = sizeof(int) * 2 + sizeof(double) * (l_max + 1) * (m_max + 1);
        in.seekg(offset_S_start, std::ios::beg);
        in.read(reinterpret_cast<char *>(S_buffer.data()), sizeof(double) * block_rows * full_row);
        for (int i = 0; i < block_rows; ++i)
        {
            for (int j = 0; j < block_cols; ++j)
            {
                S_block[i][j] = S_buffer[i * full_row + j];
            }
        }
        in.close();
        return true;
    }
    // Reads only the C and S coefficients for a specific (l, m) from the binary file
    // Returns true if successful, false otherwise
    bool read_coefficient_binary(const std::string &filename, int l, int m, double &C_val, double &S_val)
    {
        std::ifstream in(filename, std::ios::binary);
        if (!in.is_open())
            return false;
        int l_max = 0, m_max = 0;
        in.read(reinterpret_cast<char *>(&l_max), sizeof(int));
        in.read(reinterpret_cast<char *>(&m_max), sizeof(int));
        if (l > l_max || m > m_max)
        {
            in.close();
            return false;
        }
        // Calculate offset for C[l][m]
        std::streamoff offset_C = sizeof(int) * 2 + sizeof(double) * ((l * (m_max + 1)) + m);
        in.seekg(offset_C, std::ios::beg);
        in.read(reinterpret_cast<char *>(&C_val), sizeof(double));
        // Calculate offset for S[l][m]
        std::streamoff offset_S = sizeof(int) * 2 + sizeof(double) * ((l_max + 1) * (m_max + 1)) + sizeof(double) * ((l * (m_max + 1)) + m);
        in.seekg(offset_S, std::ios::beg);
        in.read(reinterpret_cast<char *>(&S_val), sizeof(double));
        in.close();
        return true;
    }
    // Writes C and S coefficient matrices to a binary file for efficient reading
    // Format: [int l_max][int m_max][C coefficients][S coefficients] (row-major)
    bool write_coefficients_binary(const std::string &filename, const std::vector<std::vector<double>> &C, const std::vector<std::vector<double>> &S, int l_max, int m_max)
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open())
            return false;
        // Write l_max and m_max
        out.write(reinterpret_cast<const char *>(&l_max), sizeof(int));
        out.write(reinterpret_cast<const char *>(&m_max), sizeof(int));
        // Write C coefficients
        for (int l = 0; l <= l_max; ++l)
        {
            for (int m = 0; m <= m_max; ++m)
            {
                double c = (l < C.size() && m < C[l].size()) ? C[l][m] : 0.0;
                out.write(reinterpret_cast<const char *>(&c), sizeof(double));
            }
        }
        // Write S coefficients
        for (int l = 0; l <= l_max; ++l)
        {
            for (int m = 0; m <= m_max; ++m)
            {
                double s = (l < S.size() && m < S[l].size()) ? S[l][m] : 0.0;
                out.write(reinterpret_cast<const char *>(&s), sizeof(double));
            }
        }
        out.close();
        return true;
    }
    // recursive normalized legendre polynomials calculates P(l,m,phi) for all degree and order up to l and m at the geocentric latitude phi. The results are stored in the 2D array P where P[l][m] corresponds to degree l and order m
    // Using Holmes and Featherstone 2002 recursion for normalized associated legendre polynomials
    // P(0,0) = 1
    // P(1,0) = sqrt(3)*sin(phi)*P(0,0)
    // P(1,1) = sqrt(3)*cos(phi)*P(0,0)
    // P(l,l) = sqrt((2l+1)/(2l))*cos(phi)*P(l-1,l-1)
    // P(l,l-1) = sqrt(2l+1)*sin(phi)*P(l-1,l-1)
    // P(l,m) = sqrt((2l+1)/(l-m)(l+m))*[sqrt(2l-1)*sin(phi)*P(l-1,m) - sqrt((l-1+m)(l-1-m)/(2l-3))*P(l-2,m)]
    std::vector<std::vector<double>> Plm_bar(int l_max, int m_max, double phi)
    {
        // theta = colatitude
        const double theta = M_PI_2 - phi;
        const double u = std::sin(theta);
        const double t = std::cos(theta);

        if (m_max > l_max)
            m_max = l_max;
        if (l_max < 0 || m_max < 0)
            return {};

        // P[n][m]
        std::vector<std::vector<double>> P(l_max + 1, std::vector<double>(m_max + 1, 0.0));

        // sectoral seeds
        P[0][0] = 1.0;
        if (l_max >= 1 && m_max >= 1)
        {
            P[1][1] = std::sqrt(3.0) * u;
        }
        for (int l = 2; l <= l_max && l <= m_max; ++l)
        {
            P[l][l] = std::sqrt((2.0 * l + 1.0) / (2.0 * l)) * u * P[l - 1][l - 1];
        }

        // columns (fixed m), increasing degree l
        for (int m = 0; m <= m_max; ++m)
        {
            if (m == 0 && l_max >= 1)
            {
                // first off-diagonal for m=0 is covered by general case below (l=m+1)
            }
            // explicit first off-diagonal: P[m+1][m] = sqrt(2m+3)*t*P[m][m]
            if (m + 1 <= l_max)
            {
                P[m + 1][m] = std::sqrt(2.0 * m + 3.0) * t * P[m][m];
            }
            for (int l = m + 2; l <= l_max; ++l)
            {
                // fully-normalized coefficients (force floating-point!)
                const double anm = std::sqrt(((2.0 * l - 1.0) * (2.0 * l + 1.0)) / ((l - m) * (l + m)));
                const double bnm = std::sqrt(((2.0 * l + 1.0) * (l + m - 1.0) * (l - m - 1.0)) /
                                             ((l - m) * (l + m) * (2.0 * l - 3.0)));
                P[l][m] = anm * t * P[l - 1][m] - bnm * P[l - 2][m];
            }
        }

        return P;
    }

    // Recursive tangent function to compute mTan(phi) = (m-1)Tan(phi) + Tan(phi).
    // returns a vector of size m_max+1 where the index corresponds to m
    std::vector<double> recursive_tangent(int m_max, double phi)
    {
        std::vector<double> mTan(m_max + 1, 0.0);
        mTan[0] = 0.0; // 0*Tan(phi) = 0
        if (m_max >= 1)
        {
            mTan[1] = tan(phi); // Tan(1*phi) = tan(phi)
        }
        for (int m = 2; m <= m_max; ++m)
        {
            mTan[m] = mTan[m - 1] + mTan[1];
        }
        return mTan;
    }

    // Recursive trig of longitude to compute SIN(m*lambda), COS(m*lambda) for m=0 to m_max. returns a vector of size m_max+1 where the index corresponds to m
    //  SIN(m*lambda) = 2*COS(lambda)*SIN((m-1)*lambda) - SIN((m-2)*lambda)
    //  COS(m*lambda) = 2*COS(lambda)*COS((m-1)*lambda) - COS((m-2)*lambda)
    std::pair<std::vector<double>, std::vector<double>> recursive_sine_cosine(int m_max, double lambda)
    {
        std::vector<double> sinL(m_max + 1, 0.0);
        std::vector<double> cosL(m_max + 1, 0.0);
        sinL[0] = 0.0; // sin(0*lambda) = 0
        cosL[0] = 1.0; // cos(0*lambda) = 1

        if (m_max >= 1)
        {
            sinL[1] = sin(lambda); // sin(1*lambda) = sin(lambda
            cosL[1] = cos(lambda); // cos(1*lambda) = cos(lambda)
        }
        for (int m = 2; m <= m_max; ++m)
        {
            sinL[m] = 2 * cosL[1] * sinL[m - 1] - sinL[m - 2];
            cosL[m] = 2 * cosL[1] * cosL[m - 1] - cosL[m - 2];
        }
        return {sinL, cosL};
    }

    // This function takes a text file of normalized spherical harmonic coefficients and loads them into the 2D vectors C and S. The file should have the following format:
    // l m C(l,m) S(l,m) {RMS error C} {RMS error S}
    //     2    0   -0.484165143790815D-03    0.000000000000000D+00    0.7481239490D-11    0.0000000000D+00
    // where l is the degree, m is the order, C(l,m) is the cosine coefficient, S(l,m) is the sine coefficient, and RMS error C and RMS error S are optional and can be ignored. The function returns true if the file was loaded successfully, and false otherwise.
    // The file is stored ordered by degree then order, so all coefficients for degree 0 are listed first, followed by all coefficients for degree 1, and so on. Within each degree, the coefficients are listed in order of increasing order.
    bool read_EGM2008_coefficients_text(const std::string &filename, int l_max, int m_max, std::vector<std::vector<double>> &C, std::vector<std::vector<double>> &S)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            return false; // Failed to open file
        }
        // Resize C and S to hold coefficients up to l_max and m_max
        C.resize(l_max + 1, std::vector<double>(m_max + 1, 0.0));
        S.resize(l_max + 1, std::vector<double>(m_max + 1, 0.0));
        std::string line;
        while (std::getline(file, line))
        {
            std::cout << "Reading line: " << line << std::endl;
            // Replace 'D' or 'd' with 'E' for scientific notation compatibility
            for (char &ch : line)
            {
                if (ch == 'D' || ch == 'd')
                    ch = 'E';
            }
            std::istringstream iss(line);
            int l, m;
            double c, s;
            if (!(iss >> l >> m >> c >> s))
            {
                continue; // Skip lines that don't match the expected format
            }
            if (l <= l_max && m <= m_max)
            {
                C[l][m] = c;
                S[l][m] = s;
            }
        }
        file.close();
        return true;
    }

    // Computes the gravitational potential at (r, phi, lambda) using spherical harmonics
    double compute_gravitational_potential(double r, double phi, double lambda, int l_max, int m_max,
                                           const std::vector<std::vector<double>> &C,
                                           const std::vector<std::vector<double>> &S,
                                           double a, double GM)
    {
        // Compute normalized Legendre polynomials
        std::vector<std::vector<double>> P = Plm_bar(l_max, m_max, phi);
        // Compute recursive sine and cosine values
        auto [sinL, cosL] = recursive_sine_cosine(m_max, lambda);
        double V = GM / r;
        for (int l = 2; l <= l_max; ++l)
        {
            double r_ratio_l = pow(a / r, l);
            for (int m = 0; m <= l; ++m)
            {
                double C_lm = C[l][m];
                double S_lm = S[l][m];
                double P_lm = P[l][m];
                double cos_m_lambda = cosL[m];
                double sin_m_lambda = sinL[m];
                V += GM / r * r_ratio_l * P_lm * (C_lm * cos_m_lambda + S_lm * sin_m_lambda);
            }
        }
        return V;
    }

    // Returns {radius, geocentric latitude (radians), longitude (radians)}
    std::array<double, 3> cartesian_to_geocentric(double x, double y, double z)
    {
        double r = std::sqrt(x * x + y * y + z * z);
        double phi_gc = std::asin(z / r); // geocentric latitude
        double lambda = std::atan2(y, x); // longitude
        if (lambda < 0)
        {
            lambda += 2 * M_PI; // Normalize to [0, 2π]
        }
        else if (lambda >= 2 * M_PI)
        {
            lambda -= 2 * M_PI; // Normalize to [0, 2π]
        }
        return {r, phi_gc, lambda};
    }
}