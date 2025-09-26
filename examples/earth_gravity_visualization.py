#!/usr/bin/env python3
"""
Earth Gravity Field Visualization using EGM2008

This script computes gravitational acceleration magnitude at evenly distributed
points on Earth's surface, accounting for Earth's flattening, and visualizes
the results using cartopy for both low and high fidelity calculations.

Requirements:
- pyshg (Spherical Harmonic Gravity library)
- numpy
- matplotlib
- cartopy
"""

import numpy as np
import matplotlib.pyplot as plt
import pyshg

# Using basic matplotlib plotting for all visualizations

def geodetic_to_geocentric(lat_geodetic_deg, lon_deg, alt_m=0.0, a=None, f=None):
    """
    Convert geodetic coordinates to geocentric coordinates.
    
    Parameters:
    lat_geodetic_deg : float or array
        Geodetic latitude in degrees
    lon_deg : float or array
        Longitude in degrees (same for both systems)
    alt_m : float or array
        Altitude above ellipsoid in meters
    a : float, optional
        Semi-major axis (default: EGM2008_A)
    f : float, optional  
        Flattening (default: WGS84 value)
        
    Returns:
    tuple : (r_geocentric, lat_geocentric_deg, lon_deg) 
        - r_geocentric: geocentric radius in meters
        - lat_geocentric_deg: geocentric latitude in degrees  
        - lon_deg: longitude in degrees (unchanged)
    """
    if a is None:
        a = pyshg.EGM2008_A
    if f is None:
        f = 1/298.257223563  # WGS84 flattening
    
    lat_geodetic_rad = np.radians(lat_geodetic_deg)
    b = a * (1 - f)  # Semi-minor axis
    
    # First eccentricity squared
    e2 = f * (2 - f)
    
    # Prime vertical radius of curvature
    N = a / np.sqrt(1 - e2 * np.sin(lat_geodetic_rad)**2)
    
    # Cartesian coordinates from geodetic
    cos_lat = np.cos(lat_geodetic_rad)
    sin_lat = np.sin(lat_geodetic_rad)
    cos_lon = np.cos(np.radians(lon_deg))
    sin_lon = np.sin(np.radians(lon_deg))
    
    X = (N + alt_m) * cos_lat * cos_lon
    Y = (N + alt_m) * cos_lat * sin_lon
    Z = (N * (1 - e2) + alt_m) * sin_lat
    
    # Convert to geocentric coordinates
    r_geocentric = np.sqrt(X**2 + Y**2 + Z**2)
    lat_geocentric_rad = np.arctan2(Z, np.sqrt(X**2 + Y**2))
    lat_geocentric_deg = np.degrees(lat_geocentric_rad)
    
    return r_geocentric, lat_geocentric_deg, lon_deg

def compute_WGS84_ellipse_gravity(phi, alt_m=0.0):
    """
    Parameters:
    phi : float or array
        Geodetic latitude in radians
    alt_m : float or array
        Altitude above ellipsoid in meters
        
    Returns:
    float or array : Normal gravitational acceleration in m/s² (with centrifugal effects)
    """
    GM = 3986004.418e8  # in m^3/s^2
    a = 6378137.0  # Semi-major axis in meters
    b = a * (1 - 1/298.257223563)  # Semi-minor axis
    f = 1/298.257223563  # Flattening
    r = np.sqrt(((a**2*np.cos(phi))**2+(b**2*np.sin(phi))**2)/((a*np.cos(phi))**2+(b*np.sin(phi))**2))
    J2 = 1.08263e-3  # Earth's second dynamic form factor
    g_normal = GM/r**2*(1-3*J2*a**2/(2*r**2)*(3*np.sin(phi)**2 - 1))
    return g_normal


def generate_grid_points(lat_resolution=2.0, lon_resolution=2.0):
    """
    Generate evenly distributed geodetic lat/lon grid points.
    
    Note: Grid points are generated in geodetic coordinates (WGS84-like),
    which is the standard for mapping and visualization. These will be
    converted to geocentric coordinates for gravity calculations.
    
    Use sphere point picking methods for more uniform distribution
    
    Parameters:
    lat_resolution : float
        Geodetic latitude spacing in degrees
    lon_resolution : float  
        Longitude spacing in degrees
        
    Returns:
    tuple : (lat_grid, lon_grid) in geodetic degrees for plotting
    """
    u = np.linspace(0, 1, int(180/lat_resolution)+1)
    v = np.linspace(0, 1, int(360/lon_resolution)+1)
    lats = np.degrees(np.arcsin(1 - 2*u))
    lons = 360 * v - 180
    
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    
    return lat_grid, lon_grid

def compute_gravity_field(lat_grid, lon_grid, max_degree, altitude_km=0.0):
    """
    Compute gravitational acceleration anomalies over a geodetic grid.
    
    This function computes the EGM2008 gravity field and subtracts the
    normal gravity field to show gravity anomalies, which reveal mass
    distribution features rather than the predictable latitude variation.
    
    Parameters:
    lat_grid : array
        Geodetic latitude grid in degrees
    lon_grid : array  
        Longitude grid in degrees (same for geodetic/geocentric)
    max_degree : int
        Maximum EGM2008 degree/order
    altitude_km : float
        Altitude above WGS84 ellipsoid in km
        
    Returns:
    array : Gravity anomalies (EGM2008 - Normal) in mGal (milliGal)
    """
    print(f"Computing gravity anomalies with degree {max_degree}...")
    print(f"Grid size: {lat_grid.shape[0]} x {lat_grid.shape[1]} = {lat_grid.size} points")
    
    gravity_anomalies = np.zeros_like(lat_grid)
    
    # Flatten for easier iteration
    lat_flat = lat_grid.flatten()
    lon_flat = lon_grid.flatten()
    
    total_points = len(lat_flat)
    
    for i, (lat_geodetic, lon) in enumerate(zip(lat_flat, lon_flat)):
        if i % 500 == 0:  # Progress indicator
            print(f"  Progress: {i}/{total_points} ({100*i/total_points:.1f}%)")
        
        # Convert geodetic coordinates to geocentric coordinates
        # lat_geodetic is in degrees, altitude_km converted to meters
        altitude_m = altitude_km * 1000
        r, lat_geocentric_deg, lon_out = geodetic_to_geocentric(
            lat_geodetic, lon, altitude_m
        )
        
        # Convert geocentric coordinates to radians for pyshg
        lat_geocentric_rad = np.radians(lat_geocentric_deg)
        lon_rad = np.radians(lon)
        lat_geodetic_rad = np.radians(lat_geodetic)
        try:
            # Compute EGM2008 gravitational acceleration vector
            g_vec = pyshg.g_EGM2008(r, lat_geocentric_rad, lon_rad, max_degree)
            g_egm2008 = np.sqrt(g_vec[0]**2 + g_vec[1]**2 + g_vec[2]**2)
            
            # Compute normal gravity at this location
            g_normal = compute_WGS84_ellipse_gravity(lat_geodetic_rad, altitude_m)
            
            # Compute gravity anomaly in milliGal (1 mGal = 10⁻⁵ m/s²)
            gravity_anomalies.flat[i] = (g_egm2008 - g_normal) * 1e5  # Convert to mGal
            
        except Exception as e:
            print(f"  Error at geodetic lat={lat_geodetic:.1f}, lon={lon:.1f}: {e}")
            gravity_anomalies.flat[i] = np.nan
    
    print(f"  Completed! Anomaly range: {np.nanmin(gravity_anomalies):.1f} to {np.nanmax(gravity_anomalies):.1f} mGal")
    
    return gravity_anomalies

def save_gravity_data_csv(lat_grid, lon_grid, gravity_data, filename):
    """
    Save gravity field data as CSV file with lat, lon, gravity columns.
    """
    try:
        import pandas as pd
        
        # Flatten the grids
        lat_flat = lat_grid.flatten()
        lon_flat = lon_grid.flatten()
        gravity_flat = gravity_data.flatten()
        
        # Create DataFrame
        df = pd.DataFrame({
            'latitude': lat_flat,
            'longitude': lon_flat,
            'gravity_anomaly_mGal': gravity_flat
        })
        
        # Remove NaN values
        df = df.dropna()
        
        # Save to CSV
        df.to_csv(filename, index=False, float_format='%.6f')
        
    except ImportError:
        # Fallback: save using numpy if pandas not available
        lat_flat = lat_grid.flatten()
        lon_flat = lon_grid.flatten()  
        gravity_flat = gravity_data.flatten()
        
        # Remove NaN entries
        valid_mask = ~np.isnan(gravity_flat)
        lat_valid = lat_flat[valid_mask]
        lon_valid = lon_flat[valid_mask]
        gravity_valid = gravity_flat[valid_mask]
        
        # Save using numpy
        header = 'latitude,longitude,gravity_anomaly_mGal'
        data = np.column_stack([lat_valid, lon_valid, gravity_valid])
        np.savetxt(filename, data, delimiter=',', header=header, 
                   comments='', fmt='%.6f')

def load_gravity_data(data_filename):
    """
    Load gravity field data from saved .npz file.
    
    Parameters:
    data_filename : str
        Path to the .npz data file
        
    Returns:
    dict : Dictionary containing lat_grid, lon_grid, gravity_anomalies, title
    """
    data = np.load(data_filename)
    return {
        'lat_grid': data['lat_grid'],
        'lon_grid': data['lon_grid'], 
        'gravity_anomalies': data['gravity_anomalies'],
        'title': str(data['title']) if 'title' in data else "Loaded Data"
    }

def plot_gravity_field_basic(lat_grid, lon_grid, gravity_mag, title, filename=None):
    """
    Plot gravity field using basic matplotlib with optimized colorbar scaling.
    """
    fig, ax = plt.subplots(figsize=(15, 10))
    
    # Use preferred colorbar range (-70 to 90 mGal) for optimal visualization
    vmin, vmax = -70, 90
    
    # Plot gravity data with smooth contours
    levels = np.linspace(vmin, vmax, 50)
    
    im = ax.contourf(lon_grid, lat_grid, gravity_mag, 
                     levels=levels, cmap='plasma', extend='both')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label='Gravity Anomalies (mGal)', shrink=0.8)
    
    # Labels and grid
    ax.set_xlabel('Longitude (degrees)', fontsize=12)
    ax.set_ylabel('Latitude (degrees)', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    
    # Add degree markings
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    
    # Data statistics
    data_min = np.nanmin(gravity_mag)
    data_max = np.nanmax(gravity_mag)
    data_mean = np.nanmean(gravity_mag)
    data_std = np.nanstd(gravity_mag)
    
    # Statistics text box
    stats_text = f"""Data Statistics:
Full Range: {data_min:.1f} to {data_max:.1f} mGal
Mean: {data_mean:.1f} ± {data_std:.1f} mGal
Display Range: {vmin} to {vmax} mGal
Coverage: {np.nanmax(lat_grid) - np.nanmin(lat_grid):.1f}° × {np.nanmax(lon_grid) - np.nanmin(lon_grid):.1f}°
Points: {gravity_mag.size:,}"""
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props, fontfamily='monospace')
    
    plt.title(f"{title}\nOptimized Range: {vmin} to {vmax} mGal", fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if filename:
        # Save the plot
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {filename}")
        
        # Save the raw data
        data_filename = filename.replace('.png', '_data.npz')
        np.savez_compressed(data_filename,
                           lat_grid=lat_grid,
                           lon_grid=lon_grid,
                           gravity_anomalies=gravity_mag,
                           title=title)
        print(f"Saved data to {data_filename}")
        
        # Save as CSV for easy access
        csv_filename = filename.replace('.png', '_data.csv')
        save_gravity_data_csv(lat_grid, lon_grid, gravity_mag, csv_filename)
        print(f"Saved CSV data to {csv_filename}")
        
    plt.show()

def analyze_gravity_differences(gravity_low, gravity_high, lat_grid, lon_grid):
    """
    Analyze differences between low and high fidelity calculations.
    """
    diff = gravity_high - gravity_low
    
    print("\n=== Gravity Anomaly Analysis ===")
    print(f"Low fidelity  - Min: {np.nanmin(gravity_low):.1f}, Max: {np.nanmax(gravity_low):.1f} mGal")
    print(f"High fidelity - Min: {np.nanmin(gravity_high):.1f}, Max: {np.nanmax(gravity_high):.1f} mGal")
    print(f"Difference    - Min: {np.nanmin(diff):.1f}, Max: {np.nanmax(diff):.1f} mGal")
    print(f"RMS difference: {np.sqrt(np.nanmean(diff**2)):.1f} mGal")
    print(f"Mean abs diff:  {np.nanmean(np.abs(diff)):.6f} m/s²")
    
    # Plot difference field using basic matplotlib
    fig, ax = plt.subplots(figsize=(15, 10))
    
    # Use symmetric colorbar around zero
    vmax = max(abs(np.nanmin(diff)), abs(np.nanmax(diff)))
    levels = np.linspace(-vmax, vmax, 50)
    
    im = ax.contourf(lon_grid, lat_grid, diff, 
                     levels=levels, cmap='RdBu_r', extend='both')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label='Gravity Difference: High - Low Fidelity (mGal)', shrink=0.8)
    
    # Labels and grid
    ax.set_xlabel('Longitude (degrees)', fontsize=12)
    ax.set_ylabel('Latitude (degrees)', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    
    plt.title('EGM2008 Gravity Field: High Fidelity - Low Fidelity Difference', fontsize=16)
    plt.tight_layout()
    plt.savefig('examples/gravity_difference_field.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """
    Main execution function.
    """
    print("=== Earth Gravity Field Visualization using EGM2008 ===")
    print()
    
    # Check if pyshg is available with EGM2008 functions
    if not hasattr(pyshg, 'g_EGM2008'):
        print("Error: pyshg module does not have g_EGM2008 function.")
        print("Please ensure you have the latest version of the pyshg module.")
        return
    
    print(f"Using EGM2008 constants:")
    print(f"  GM = {pyshg.EGM2008_GM:,.0f} m³/s²")
    print(f"  A  = {pyshg.EGM2008_A:,.1f} m")
    print()
    
    # Configuration
    configs = [
        # {
        #     'name': 'Low Fidelity',
        #     'degree': 36,  # Degree 36 ≈ 550 km resolution
        #     'lat_res': 10.0,  # 10 degree grid for faster computation
        #     'lon_res': 10.0,
        #     'filename': 'examples/gravity_anomalies_low_fidelity.png'
        # },
        {
            'name': 'High Fidelity', 
            'degree': 1800,  # Degree 120 ≈ 165 km resolution
            'lat_res': 1,  # 5 degree grid  
            'lon_res': 1,
            'filename': 'examples/gravity_anomalies_high_fidelity.png'
        }
    ]
    
    results = {}
    
    for config in configs:
        print(f"\n=== {config['name']} Calculation (Degree {config['degree']}) ===")
        
        # Generate grid
        lat_grid, lon_grid = generate_grid_points(config['lat_res'], config['lon_res'])
        
        # Compute gravity anomalies
        gravity_anomalies = compute_gravity_field(lat_grid, lon_grid, config['degree'])
        
        # Plot
        title = f"Earth Gravity Anomalies - {config['name']} (EGM2008 Degree {config['degree']})"
        plot_gravity_field_basic(lat_grid, lon_grid, gravity_anomalies, title, config['filename'])
        
        # Store results
        results[config['name']] = {
            'gravity': gravity_anomalies,
            'lat_grid': lat_grid,
            'lon_grid': lon_grid,
            'degree': config['degree']
        }
    
    # Analyze differences if we have both results
    if len(results) == 2:
        # Interpolate low fidelity to high fidelity grid for comparison
        from scipy.interpolate import griddata
        
        low_res = results['Low Fidelity']
        high_res = results['High Fidelity'] 
        
        # Flatten low resolution data
        points_low = np.column_stack([low_res['lat_grid'].flatten(), 
                                      low_res['lon_grid'].flatten()])
        values_low = low_res['gravity'].flatten()
        
        # Interpolate to high resolution grid
        points_high = np.column_stack([high_res['lat_grid'].flatten(),
                                       high_res['lon_grid'].flatten()])
        
        gravity_low_interp = griddata(points_low, values_low, points_high, method='linear')
        gravity_low_interp = gravity_low_interp.reshape(high_res['gravity'].shape)
        
        # analyze_gravity_differences(gravity_low_interp, high_res['gravity'], 
                                    # high_res['lat_grid'], high_res['lon_grid'])
    
    print("\n=== Gravity Anomaly Visualization Complete ===")
    print("Generated files:")
    for config in configs:
        print(f"  {config['filename']}")
    if len(results) == 2:
        print("  examples/gravity_anomaly_difference_field.png")
    
    # Send completion notification
    try:
        from examples.notification_utils import notify_completion
        file_list = [config['filename'] for config in configs]
        notify_completion(f"EGM2008 gravity visualization complete! Generated {len(file_list)} files.")
    except ImportError:
        print("Notification utility not found; skipping notification.")

if __name__ == "__main__":
    main()