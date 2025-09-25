/*
 * Python bindings for Spherical Harmonic Gravity (SHG) library
 * -------------------------------------------------------------
 * 
 * This file provides Python bindings for the SHG C++ library using pybind11.
 * It exposes the main gravitational acceleration and potential functions to Python.
 *
 * Project: Spherical-Harmonic-Gravity
 * License: MIT (see LICENSE file for details)
 * Author: DavidMSCode
 * Created: 2025
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <array>
#include "SHG.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyshg, m) {
    m.doc() = "Python bindings for Spherical Harmonic Gravity (SHG) library";

    // Wrapper functions for template functions that return std::array<double, 3>
    auto compute_acceleration_array = [](double r, double phi, double lambda, int l_max, int m_max,
                                        const std::vector<std::vector<double>>& C,
                                        const std::vector<std::vector<double>>& S,
                                        double a, double GM) -> std::array<double, 3> {
        return SHG::g<std::array<double, 3>>(r, phi, lambda, l_max, m_max, C, S, a, GM);
    };

    auto compute_acceleration_from_cartesian_array = [](double x, double y, double z, int l_max, int m_max,
                                                       const std::vector<std::vector<double>>& C,
                                                       const std::vector<std::vector<double>>& S,
                                                       double a, double GM) -> std::array<double, 3> {
        return SHG::gravitational_acceleration_from_cartesian<std::array<double, 3>>(x, y, z, l_max, m_max, C, S, a, GM);
    };

    // Main functions for gravitational computation
    m.def("g", compute_acceleration_array,
          "Compute gravitational acceleration at geocentric coordinates (r, phi, lambda)",
          py::arg("r"), py::arg("phi"), py::arg("lambda"), py::arg("l_max"), py::arg("m_max"),
          py::arg("C"), py::arg("S"), py::arg("a"), py::arg("GM"));

    m.def("gravitational_acceleration_from_cartesian", compute_acceleration_from_cartesian_array,
          "Compute gravitational acceleration from Cartesian coordinates (x, y, z)",
          py::arg("x"), py::arg("y"), py::arg("z"), py::arg("l_max"), py::arg("m_max"),
          py::arg("C"), py::arg("S"), py::arg("a"), py::arg("GM"));

    m.def("pot", &SHG::pot,
          "Compute gravitational potential at geocentric coordinates (r, phi, lambda)",
          py::arg("r"), py::arg("phi"), py::arg("lambda"), py::arg("l_max"), py::arg("m_max"),
          py::arg("C"), py::arg("S"), py::arg("a"), py::arg("GM"));

    // Utility functions
    m.def("cartesian_to_geocentric", &SHG::cartesian_to_geocentric,
          "Convert Cartesian coordinates to geocentric (r, phi, lambda)",
          py::arg("x"), py::arg("y"), py::arg("z"));

    m.def("Plm_bar", &SHG::Plm_bar,
          "Compute normalized associated Legendre polynomials",
          py::arg("l_max"), py::arg("m_max"), py::arg("phi"));

    m.def("recursive_tangent", &SHG::recursive_tangent,
          "Compute vector of m*Tan(phi) for m=0 to m_max",
          py::arg("m_max"), py::arg("phi"));

    m.def("recursive_sine_cosine", &SHG::recursive_sine_cosine,
          "Compute vectors of sin(m*lambda) and cos(m*lambda) for m=0 to m_max",
          py::arg("m_max"), py::arg("lambda"));

    // Coefficient I/O functions - wrapper that returns the matrices
    auto read_coefficients_wrapper = [](const std::string& filename, int l, int m) -> py::tuple {
        std::vector<std::vector<double>> C_block, S_block;
        bool success = SHG::read_coefficients_block_binary(filename, l, m, C_block, S_block);
        return py::make_tuple(success, C_block, S_block);
    };
    
    m.def("read_coefficients_block_binary_new", read_coefficients_wrapper,
          "Read spherical harmonic coefficients from binary file - returns (success, C, S)",
          py::arg("filename"), py::arg("l"), py::arg("m"));
          
    // Keep the old function for backward compatibility (but it doesn't work properly)
    m.def("read_coefficients_block_binary", &SHG::read_coefficients_block_binary,
          "Read spherical harmonic coefficients from binary file (deprecated - use read_coefficients_block_binary_new)",
          py::arg("filename"), py::arg("l"), py::arg("m"), py::arg("C_block"), py::arg("S_block"));

    m.def("read_coefficient_binary", &SHG::read_coefficient_binary,
          "Read a specific spherical harmonic coefficient from binary file",
          py::arg("filename"), py::arg("l"), py::arg("m"), py::arg("C_val"), py::arg("S_val"));

    m.def("write_coefficients_binary", &SHG::write_coefficients_binary,
          "Write spherical harmonic coefficients to binary file",
          py::arg("filename"), py::arg("C"), py::arg("S"), py::arg("l_max"), py::arg("m_max"));

    m.def("read_EGM2008_coefficients_text", &SHG::read_EGM2008_coefficients_text,
          "Read EGM2008 spherical harmonic coefficients from text file",
          py::arg("filename"), py::arg("l_max"), py::arg("m_max"), py::arg("C"), py::arg("S"));

    // EGM2008-specific functions
    m.def("g_EGM2008", &SHG::g_EGM2008,
          "Compute gravitational acceleration using EGM2008 model with user-specified degree",
          py::arg("r"), py::arg("phi"), py::arg("lambda"), 
          py::arg("max_degree"), py::arg("coefficient_path") = "");
          
    m.def("U_EGM2008", &SHG::U_EGM2008,
          "Compute gravitational potential using EGM2008 model with user-specified degree",
          py::arg("r"), py::arg("phi"), py::arg("lambda"), 
          py::arg("max_degree"), py::arg("coefficient_path") = "");
          
    m.def("load_EGM2008_coefficients", &SHG::load_EGM2008_coefficients,
          "Load EGM2008 coefficients with automatic fallback from binary to text format",
          py::arg("coefficient_path"), py::arg("C"), py::arg("S"));

    // Constants from SHG header
    m.attr("WGS84_A") = SHG::WGS84_A;           // WGS84 semi-major axis in meters
    m.attr("WGS84_A_KM") = SHG::WGS84_A_KM;     // WGS84 semi-major axis in kilometers
    m.attr("WGS84_GM") = SHG::WGS84_GM;         // WGS84 gravitational parameter in m³/s²
    m.attr("WGS84_GM_KM") = SHG::WGS84_GM_KM;   // WGS84 gravitational parameter in km³/s²
    m.attr("WGS84_F") = SHG::WGS84_F;           // WGS84 flattening
    
    m.attr("EGM2008_A") = SHG::EGM2008_A;       // EGM2008 semi-major axis in meters
    m.attr("EGM2008_A_KM") = SHG::EGM2008_A_KM; // EGM2008 semi-major axis in kilometers
    m.attr("EGM2008_GM") = SHG::EGM2008_GM;     // EGM2008 gravitational parameter in m³/s²
    m.attr("EGM2008_GM_KM") = SHG::EGM2008_GM_KM; // EGM2008 gravitational parameter in km³/s²
    m.attr("EGM2008_MAX_DEGREE") = SHG::EGM2008_MAX_DEGREE; // EGM2008 maximum degree
    m.attr("EGM2008_MAX_ORDER") = SHG::EGM2008_MAX_ORDER;   // EGM2008 maximum order
}