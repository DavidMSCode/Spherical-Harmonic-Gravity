#!/usr/bin/env python3
"""
Example usage of the Python SHG (Spherical Harmonic Gravity) library

This example demonstrates how to:
1. Load spherical harmonic coefficients
2. Compute gravitational acceleration and potential
3. Compare results from different coordinate systems
"""

import numpy as np
import pyshg

def example_multiple_points():
    """Example with multiple test points"""
    print("\n=== Multiple Points Example ===")
    
    # Load coefficients using the new working function
    l_max, m_max = 100, 100
    
    # Use EGM2008 constants - EGM2008 functions expect meters, not km
    a = pyshg.EGM2008_A     # EGM2008 reference radius in meters
    
    # Test points at different latitudes and altitudes
    altitudes = [0]  # km (will be converted to meters)
    latitudes = [-90,-45,0.0,45,90]       # degrees
    
    print(f"{'Radius (km)':>18} {'Altitude (km)':>12} {'Latitude (°)':>12} {'Accel Mag (m/s²)':>18} {'Potential (m²/s²)':>18}")
    print("-" * 90)
    
    for alt in altitudes:
        for lat in latitudes:
            # Compute geocentric radius accounting for flattening (WGS84)
            f = 1/298.257223563
            b = a * (1 - f)
            phi = np.radians(lat)
            #get surface radius at latitude
            cos_phi = np.cos(phi)
            sin_phi = np.sin(phi)
            r = np.sqrt((a**2*cos_phi**2 + b**2*sin_phi**2)) + alt * 1000  # Convert km to meters
            phi = np.radians(lat)
            lambda_ = 0.0
            
            # EGM2008 functions expect meters and radians, return m/s² and m²/s²
            accel = pyshg.g_EGM2008(r, phi, lambda_, l_max)
            potential = pyshg.U_EGM2008(r, phi, lambda_, l_max)
            
            accel_mag = np.linalg.norm(accel)  # Already in m/s²
            
            print(f"{r/1000:18.9f}{alt:12.0f} {lat:12.0f} {accel_mag:18.9f} {potential:18.3f}")

if __name__ == "__main__":
    example_multiple_points()
    print("\n=== SHG Python Examples Complete ===")