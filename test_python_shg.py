#!/usr/bin/env python3
"""
Test script for the Python SHG (Spherical Harmonic Gravity) library

This script tests the EGM2008-specific functions that automatically handle
coefficient loading and use the same interface as the C++ calculators.
"""

import numpy as np
import pyshg

def test_egm2008_functions():
    """Test the EGM2008 functions with various coordinates"""
    print("=== Testing EGM2008 Python Functions ===")
    print(f"EGM2008 Constants:")
    print(f"  GM = {pyshg.EGM2008_GM:,.0f} m³/s²")
    print(f"  A  = {pyshg.EGM2008_A:,.1f} m") 
    print(f"  Max Degree = {pyshg.EGM2008_MAX_DEGREE}")
    print()
    
    # Test points: same as used in C++ calculator tests
    test_points = [
        (6378137, 0.0, 0.0, 360, "Equator, sea level"),
        (6378237, 38.9072, -77.0369, 180, "Washington DC area"),
        (6378137, 45.0, 90.0, 180, "45°N, 90°E"),
        (6378137, -45.0, 0.0, 180, "45°S, 0°E"),
        (6378137, 90.0, 0.0, 180, "North pole"),
    ]
    
    print(f"{'Location':<25} {'r (m)':<10} {'φ (°)':<8} {'λ (°)':<10} {'Degree':<7} {'Potential (m²/s²)':<18} {'Gravity (m/s²)':<15}")
    print("-" * 110)
    
    for r, phi, lam, degree, description in test_points:
        # Convert to radians for computation
        phi_rad = np.radians(phi)
        lam_rad = np.radians(lam)
        
        # Compute potential and gravity
        potential = pyshg.U_EGM2008(r, phi_rad, lam_rad, degree)
        gravity_vec = pyshg.g_EGM2008(r, phi_rad, lam_rad, degree)
        gravity_mag = np.sqrt(gravity_vec[0]**2 + gravity_vec[1]**2 + gravity_vec[2]**2)
        
        print(f"{description:<25} {r:<10} {phi:<8.3f} {lam:<10.3f} {degree:<7} {potential:<18.6f} {gravity_mag:<15.6f}")

def test_coordinate_conversions():
    """Test coordinate conversion functions"""
    print("\n=== Testing Coordinate Conversions ===")
    
    # Test point
    x, y, z = 4000000.0, 3000000.0, 5000000.0  # Cartesian coordinates in meters
    
    print(f"Cartesian coordinates: x={x/1000:.1f} km, y={y/1000:.1f} km, z={z/1000:.1f} km")
    
    # Convert to geocentric
    r, phi, lam = pyshg.cartesian_to_geocentric(x, y, z)
    
    print(f"Geocentric coordinates: r={r/1000:.1f} km, φ={np.degrees(phi):.3f}°, λ={np.degrees(lam):.3f}°")
    
    # Convert back to verify
    x_back = r * np.cos(phi) * np.cos(lam)
    y_back = r * np.cos(phi) * np.sin(lam)
    z_back = r * np.sin(phi)
    
    print(f"Converted back: x={x_back/1000:.1f} km, y={y_back/1000:.1f} km, z={z_back/1000:.1f} km")
    
    # Check differences
    dx = abs(x - x_back)
    dy = abs(y - y_back) 
    dz = abs(z - z_back)
    
    print(f"Conversion errors: dx={dx:.2e} m, dy={dy:.2e} m, dz={dz:.2e} m")

def test_gravity_from_cartesian():
    """Test gravitational acceleration from Cartesian coordinates"""
    print("\n=== Testing Gravity from Cartesian Coordinates ===")
    
    # Test point above Washington DC
    lat_deg = 38.9072
    lon_deg = -77.0369
    alt = 100000  # 100 km altitude
    
    # Convert to Cartesian (approximate)
    lat_rad = np.radians(lat_deg)
    lon_rad = np.radians(lon_deg)
    r = pyshg.EGM2008_A + alt
    
    x = r * np.cos(lat_rad) * np.cos(lon_rad)
    y = r * np.cos(lat_rad) * np.sin(lon_rad)
    z = r * np.sin(lat_rad)
    
    print(f"Test location: {lat_deg:.3f}°N, {lon_deg:.3f}°E, {alt/1000:.0f} km altitude")
    print(f"Cartesian: x={x/1000:.1f} km, y={y/1000:.1f} km, z={z/1000:.1f} km")
    
    # Method 1: Direct geocentric calculation
    r_geo, phi_geo, lam_geo = pyshg.cartesian_to_geocentric(x, y, z)
    gravity1 = pyshg.g_EGM2008(r_geo, phi_geo, lam_geo, 360)
    
    # Method 2: Use the generic gravity function with loaded coefficients 
    try:
        # Load a small set of coefficients for comparison
        l_max, m_max = 10, 10
        success, C, S = pyshg.read_coefficients_block_binary_new("EGM2008Coeffs.bin", l_max, m_max)
        if success:
            gravity2 = pyshg.gravitational_acceleration_from_cartesian(
                x/1000, y/1000, z/1000,  # Convert to km for generic function
                l_max, m_max, C, S, 
                pyshg.EGM2008_A_KM, pyshg.EGM2008_GM_KM
            )
            # Convert back to m/s²
            gravity2 = [g * 1000 for g in gravity2]
            
            print(f"EGM2008 gravity (degree 360): [{gravity1[0]:.6f}, {gravity1[1]:.6f}, {gravity1[2]:.6f}] m/s²")
            print(f"Generic gravity (degree {l_max}):  [{gravity2[0]:.6f}, {gravity2[1]:.6f}, {gravity2[2]:.6f}] m/s²")
            
            mag1 = np.sqrt(sum(g**2 for g in gravity1))
            mag2 = np.sqrt(sum(g**2 for g in gravity2))
            
            print(f"Magnitude comparison: EGM2008={mag1:.6f} m/s², Generic={mag2:.6f} m/s²")
            print(f"Difference: {abs(mag1-mag2):.6f} m/s² ({abs(mag1-mag2)/mag1*100:.3f}%)")
        else:
            print("Could not load coefficients for comparison")
            
    except Exception as e:
        print(f"Coefficient loading test failed: {e}")
        print(f"EGM2008 gravity only: [{gravity1[0]:.6f}, {gravity1[1]:.6f}, {gravity1[2]:.6f}] m/s²")

def compare_with_cpp_calculators():
    """Compare Python results with C++ calculator results"""
    print("\n=== Comparing with C++ Calculator Results ===")
    
    # Use same test points as C++ tests
    test_cases = [
        (6378137, 0, 0, 360),
        (6378137, 45, 90, 180)
    ]
    
    for r, phi_deg, lam_deg, degree in test_cases:
        phi_rad = np.radians(phi_deg)
        lam_rad = np.radians(lam_deg)
        
        potential = pyshg.U_EGM2008(r, phi_rad, lam_rad, degree)
        gravity = pyshg.g_EGM2008(r, phi_rad, lam_rad, degree)
        gravity_mag = np.sqrt(sum(g**2 for g in gravity))
        
        print(f"r={r}m, φ={phi_deg}°, λ={lam_deg}°, degree={degree}")
        print(f"  Potential: {potential:.6f} m²/s²")
        print(f"  Gravity:   {gravity_mag:.6f} m/s²")
        print(f"  Components: [{gravity[0]:.6f}, {gravity[1]:.6f}, {gravity[2]:.6f}] m/s²")
        print()

if __name__ == "__main__":
    test_egm2008_functions()
    test_coordinate_conversions() 
    test_gravity_from_cartesian()
    compare_with_cpp_calculators()
    print("=== Python SHG Tests Complete ===")