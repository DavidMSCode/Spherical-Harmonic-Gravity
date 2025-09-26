#!/usr/bin/env python3
"""
Example script to load and work with saved gravity field data.

This demonstrates how to:
1. Load saved gravity data from .npz files
2. Load CSV data for use in other applications
3. Re-plot the data
4. Export to other formats
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colorbar, colors

import sys
import os
from datetime import datetime

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the load function from the visualization module
import importlib.util
spec = importlib.util.spec_from_file_location("earth_gravity_viz", 
    os.path.join(os.path.dirname(__file__), "earth_gravity_visualization.py"))
earth_gravity_viz = importlib.util.module_from_spec(spec)
spec.loader.exec_module(earth_gravity_viz)

# Make the load function available
load_gravity_data = earth_gravity_viz.load_gravity_data

# Try to import cartopy
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    CARTOPY_AVAILABLE = True
except ImportError:
    CARTOPY_AVAILABLE = False
    print("Warning: cartopy not available. Install with: conda install cartopy")

def plot_gravity_field_basic(lat_grid, lon_grid, gravity_data, title):
    """Create a simple cartopy plot or fallback to basic matplotlib"""
    
    # Use preferred colorbar range
    vmin, vmax = -70, 90
    
    if CARTOPY_AVAILABLE:
        # Simple Mollweide projection cartopy plot
        fig = plt.figure(figsize=(15, 8))
        ax = plt.axes(projection=ccrs.Mollweide())
        
        # Plot data with simple contours (transform from PlateCarree)
        im = ax.contourf(lon_grid, lat_grid, gravity_data, 
                        levels=30, vmin=vmin, vmax=vmax,
                        cmap='plasma', extend='both',
                        transform=ccrs.PlateCarree())
        
        # Add coastlines
        ax.coastlines()
        
        # Add gridlines (no labels for Mollweide - they don't work well)
        ax.gridlines(alpha=0.3)
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Gravity Anomalies (mGal)')
        
    else:
        # Basic matplotlib fallback
        fig, ax = plt.subplots(figsize=(12, 8))
        
        im = ax.contourf(lon_grid, lat_grid, gravity_data,
                        levels=30, vmin=vmin, vmax=vmax,
                        cmap='plasma', extend='both')
        
        ax.set_xlabel('Longitude (degrees)')
        ax.set_ylabel('Latitude (degrees)')
        ax.grid(True, alpha=0.3)
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Gravity Anomalies (mGal)')
    
    plt.title(f"{title}\nRange: {vmin} to {vmax} mGal")
    plt.tight_layout()
    plt.show()



def demonstrate_data_loading():
    """Demonstrate loading and using saved gravity data"""
    
    print("=== Gravity Data Loading Demo ===")
    
    # List available data files
    data_files = []
    for file in os.listdir('examples/'):
        if file.endswith('_data.npz'):
            data_files.append(os.path.join('examples/', file))
    
    if not data_files:
        print("No gravity data files found. Run the main visualization script first.")
        return
    
    print(f"Found {len(data_files)} data files:")
    for i, file in enumerate(data_files):
        print(f"  {i+1}. {file}")
    
    # Let user choose which file to load interactively
    if data_files:
        if len(data_files) > 1:
            while True:
                try:
                    choice = input(f"\nEnter file number (1-{len(data_files)}) for interactive plot [1]: ").strip()
                    if choice == "":
                        choice = "1"
                    file_idx = int(choice) - 1
                    if 0 <= file_idx < len(data_files):
                        break
                    else:
                        print(f"Please enter a number between 1 and {len(data_files)}")
                except ValueError:
                    print("Please enter a valid number")
        else:
            file_idx = 0
            
        data_file = data_files[file_idx]
        print(f"\nLoading: {data_file}")
        
        # Load the data
        data = load_gravity_data(data_file)
        
        print(f"Data shape: {data['gravity_anomalies'].shape}")
        print(f"Lat range: {np.nanmin(data['lat_grid']):.1f}° to {np.nanmax(data['lat_grid']):.1f}°")
        print(f"Lon range: {np.nanmin(data['lon_grid']):.1f}° to {np.nanmax(data['lon_grid']):.1f}°")
        print(f"Anomaly range: {np.nanmin(data['gravity_anomalies']):.1f} to {np.nanmax(data['gravity_anomalies']):.1f} mGal")
        
        # Create plot with optimized colorbar range (-70 to 90 mGal)
        plot_gravity_field_basic(
            data['lat_grid'],
            data['lon_grid'],
            data['gravity_anomalies'],
            title=f"Gravity Field: {os.path.basename(data_file)}"
        )
        
        # Save the replotted version
        replot_filename = data_file.replace('_data.npz', '_optimized_range.png')
        print(f"Optimized range plot saved to: {replot_filename}")

def export_to_other_formats():
    """Export gravity data to various formats"""
    
    print("\n=== Exporting to Other Formats ===")
    
    # Find CSV files
    csv_files = []
    for file in os.listdir('examples/'):
        if file.endswith('_data.csv'):
            csv_files.append(os.path.join('examples/', file))
    
    if not csv_files:
        print("No CSV data files found.")
        return
        
    csv_file = csv_files[0]
    print(f"Working with: {csv_file}")
    
    # Load CSV data
    try:
        import pandas as pd
        df = pd.read_csv(csv_file)
        print(f"Loaded {len(df)} data points")
        
        # Export to different formats
        base_name = csv_file.replace('_data.csv', '')
        
        # Excel format
        excel_file = f"{base_name}_data.xlsx"
        df.to_excel(excel_file, index=False)
        print(f"Exported to Excel: {excel_file}")
        
        # JSON format
        json_file = f"{base_name}_data.json"
        df.to_json(json_file, orient='records', indent=2)
        print(f"Exported to JSON: {json_file}")
        
        # Parquet format (efficient for large datasets)
        try:
            parquet_file = f"{base_name}_data.parquet"
            df.to_parquet(parquet_file, index=False)
            print(f"Exported to Parquet: {parquet_file}")
        except ImportError:
            print("Parquet export skipped (pyarrow not installed)")
            
        # Show data statistics
        print(f"\nData Statistics:")
        print(df.describe())
        
    except ImportError:
        print("pandas not available - using numpy to read CSV")
        data = np.loadtxt(csv_file, delimiter=',', skiprows=1)
        print(f"Loaded {len(data)} data points using numpy")

if __name__ == "__main__":
    demonstrate_data_loading()
    export_to_other_formats()