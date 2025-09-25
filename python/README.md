# Python Bindings for SHG (Spherical Harmonic Gravity) Library

This directory contains Python bindings for the SHG C++ library, allowing you to compute gravitational potential and acceleration using spherical harmonic coefficients directly from Python.

## Installation

### Prerequisites

- Python 3.6 or later
- C++17 compatible compiler
- CMake (for building the C++ library)

### Install Python Dependencies

```bash
pip install -r requirements.txt
```

### Build and Install the Python Module

From the root directory of the project:

```bash
# Build the Python extension
python setup.py build_ext --inplace

# Or install it system-wide
pip install .

# Or install in development mode
pip install -e .
```

## Quick Start

```python
import numpy as np
import pyshg

# Create simple test coefficients
l_max, m_max = 2, 2
C = [[0.0 for _ in range(m_max + 1)] for _ in range(l_max + 1)]
S = [[0.0 for _ in range(m_max + 1)] for _ in range(l_max + 1)]

# Set J2 coefficient
C[2][0] = -1.082635854e-3

# Earth parameters
a = pyshg.WGS84_A        # WGS84 semi-major axis
GM = pyshg.WGS84_GM      # WGS84 gravitational parameter

# Test point: 400 km altitude above equator
r = a + 400000.0         # meters
phi = 0.0               # latitude (radians)
lambda_ = 0.0           # longitude (radians)

# Compute gravitational acceleration
accel = pyshg.compute_gravitational_acceleration(r, phi, lambda_, l_max, m_max, C, S, a, GM)
print(f"Acceleration: {accel} m/s²")

# Compute gravitational potential
potential = pyshg.compute_gravitational_potential(r, phi, lambda_, l_max, m_max, C, S, a, GM)
print(f"Potential: {potential} m²/s²")
```

## Available Functions

### Main Computation Functions

- `compute_gravitational_acceleration(r, phi, lambda, l_max, m_max, C, S, a, GM)` - Compute gravitational acceleration at geocentric coordinates
- `gravitational_acceleration_from_cartesian(x, y, z, l_max, m_max, C, S, a, GM)` - Compute gravitational acceleration from Cartesian coordinates  
- `compute_gravitational_potential(r, phi, lambda, l_max, m_max, C, S, a, GM)` - Compute gravitational potential

### Utility Functions

- `cartesian_to_geocentric(x, y, z)` - Convert Cartesian to geocentric coordinates
- `Plm_bar(l_max, m_max, phi)` - Compute normalized associated Legendre polynomials
- `recursive_tangent(m_max, phi)` - Compute m*tan(phi) recursively
- `recursive_sine_cosine(m_max, lambda)` - Compute sin/cos(m*lambda) recursively

### File I/O Functions

- `read_coefficients_block_binary(filename, l, m, C_block, S_block)` - Read coefficients from binary file
- `read_coefficient_binary(filename, l, m, C_val, S_val)` - Read single coefficient from binary file
- `write_coefficients_binary(filename, C, S, l_max, m_max)` - Write coefficients to binary file
- `read_EGM2008_coefficients_text(filename, l_max, m_max, C, S)` - Read EGM2008 text format

### Constants

- `WGS84_A` - WGS84 semi-major axis (6378137.0 m)
- `WGS84_GM` - WGS84 gravitational parameter (3.986004418e14 m³/s²)
- `EGM2008_GM` - EGM2008 gravitational parameter (3.986004415e14 m³/s²)

## Function Parameters

### Coordinate Systems

- **Geocentric coordinates**: `(r, phi, lambda)` where:
  - `r` - radial distance from Earth's center (m)
  - `phi` - geocentric latitude (radians)
  - `lambda` - longitude (radians)

- **Cartesian coordinates**: `(x, y, z)` in Earth-centered, Earth-fixed (ECEF) frame (m)

### Spherical Harmonic Parameters

- `l_max` - maximum degree of spherical harmonic expansion
- `m_max` - maximum order of spherical harmonic expansion
- `C` - 2D list of cosine coefficients C[l][m]
- `S` - 2D list of sine coefficients S[l][m]
- `a` - reference radius (typically Earth's semi-major axis)
- `GM` - gravitational parameter

## Example Usage

Run the example script:

```bash
python python/example_usage.py
```

This demonstrates:
- Basic gravitational computation
- Multiple test points
- Utility function usage
- Coordinate system conversions

## Performance Notes

- The C++ implementation uses efficient recursive algorithms for Associated Legendre Functions
- For high-degree/order models (l,m > 1900), numerical stability may become an issue
- Binary coefficient files are recommended for large gravity models like EGM2008

## Troubleshooting

### Build Issues

If you encounter build issues:

1. Make sure you have a C++17 compatible compiler
2. Verify pybind11 is installed: `pip install pybind11`
3. Try building with verbose output: `python setup.py build_ext --inplace --verbose`

### Import Issues

If you get import errors after building:
1. Make sure the build was successful
2. Try installing in development mode: `pip install -e .`
3. Check that the shared library was created in the correct location

### Numerical Issues

For high-degree spherical harmonic models:
1. Use appropriate precision for your application
2. Consider the numerical stability limits of the forward column method
3. Validate results against known benchmarks (e.g., EGM2008)

## License

MIT License - see the main project LICENSE file for details.