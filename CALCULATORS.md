# EGM2008 Calculator Programs

This directory contains two calculator programs for computing gravitational quantities using the EGM2008 Earth Gravitational Model. Both programs support interactive and command-line modes.

## Programs

### 1. GravityCalc - Gravitational Acceleration Calculator
Computes the three components of gravitational acceleration vector and its magnitude.

**Interactive Mode:**
```bash
./build/GravityCalc
```

**Command-Line Mode:**
```bash
./build/GravityCalc [--verbose|-v] <r> <phi> <lambda> <degree> [coeff_path]
./build/GravityCalc --help    # Show usage information
```

**Example:**
```bash
./build/GravityCalc 6378137 0.0 0.0 360
# Output: 9.814273

# With verbose coefficient loading messages:
./build/GravityCalc --verbose 6378137 0.0 0.0 360
```

**Parameters:**
- `--verbose, -v`: Show coefficient loading messages (optional)
- `r`: Geocentric radius in meters
- `phi`: Geocentric latitude in degrees (-90 to +90)
- `lambda`: Longitude in degrees (-180 to +180)  
- `degree`: Maximum degree/order (2 to 2190, recommend <=1900)
- `coeff_path`: Optional path to EGM2008 coefficient files

**Interactive Outputs:**
- Gravitational acceleration components (gI, gJ, gK) in m/s^2
- Acceleration magnitude in m/s^2
- Percentage difference from standard gravity (9.80665 m/s^2)

**Command-Line Output:**
- Single value: gravitational acceleration magnitude in m/s^2

### 2. PotentialCalc - Gravitational Potential Calculator
Computes the gravitational potential and related quantities.

**Interactive Mode:**
```bash
./build/PotentialCalc
```

**Command-Line Mode:**
```bash
./build/PotentialCalc [--verbose|-v] <r> <phi> <lambda> <degree> [coeff_path]
./build/PotentialCalc --help    # Show usage information
```

**Example:**
```bash
./build/PotentialCalc 6378137 0.0 0.0 360
# Output: 62528864.8

# With verbose coefficient loading messages:
./build/PotentialCalc --verbose 6378137 0.0 0.0 360
```

**Parameters:**
- `--verbose, -v`: Show coefficient loading messages (optional)
- `r`: Geocentric radius in meters
- `phi`: Geocentric latitude in degrees (-90 to +90)
- `lambda`: Longitude in degrees (-180 to +180)  
- `degree`: Maximum degree/order (2 to 2190, recommend <=1900)
- `coeff_path`: Optional path to EGM2008 coefficient files

**Interactive Outputs:**
- Gravitational potential in m^2/s^2

**Command-Line Output:**
- Single value: gravitational potential in m^2/s^2

## Building

```bash
mkdir -p build
cd build
cmake ..
make
```

## EGM2008 Coefficient Files

The programs require EGM2008 coefficient files:
- **Binary format**: `EGM2008Coeffs.bin` (faster loading)
- **Text format**: `EGM2008_to2190_TideFree.txt` (automatic fallback)

Download from: https://earth-info.nga.mil/php/download.php?file=egm-08spherical

## Coordinate System

Both calculators use **geocentric coordinates**:
- `r`: Geocentric radius in meters (distance from Earth's center)
- `phi`: Geocentric latitude in degrees (-90° to +90°)
- `lambda`: Longitude in degrees (-180° to +180°)

**Note:** This is different from geodetic coordinates. For approximate conversion:
- Geocentric radius ≈ 6,378,137 m (Earth's equatorial radius) + altitude
- Geocentric latitude ≈ geodetic latitude (small difference due to Earth's oblate shape)

## Example Session

```
Enter geocentric radius (meters): 6378237
Enter geocentric latitude (degrees): 38.9072
Enter longitude (degrees): -77.0369
Enter max degree (2-2190, recommend <=1900): 180
Enter coefficient file path (or press Enter for default): 

Computing gravitational acceleration...

========== RESULTS ==========
Location: (r=6378237.000000m, φ=38.907200°, λ=-77.036900°)
EGM2008 degree/order: 180

Gravitational acceleration components:
  gI (radial):         -9.796159 m/s²
  gJ (south):           0.000241 m/s²
  gK (east):           -0.000189 m/s²

Magnitude:              9.796159 m/s²

Difference from standard gravity (9.80665 m/s^2): -0.118%
```

## Verbosity Control

Both calculators support verbosity control:
- **Interactive mode**: Always shows coefficient loading messages and detailed results
- **Command-line mode**: Silent by default (shows only computed values)
- **Verbose flag**: Use `--verbose` or `-v` to show coefficient loading messages in command-line mode

This provides clean, automation-friendly output by default while still allowing verbose diagnostics when needed.

## Notes

- **Coordinate System**: Uses geocentric coordinates where gI points radially outward, gJ points south, gK points east
- **Performance**: Lower degrees compute faster. Degree 180 is ~3x faster than 2190 with minimal accuracy loss for most applications
- **Stability**: Degrees > 1900 may have numerical instability. A warning will be displayed.
- **Accuracy**: Full EGM2008 (degree 2190) provides centimeter-level gravitational accuracy