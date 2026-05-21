# Spherical-Harmonic-Gravity

Compute gravitational acceleration and potential from arbitrary spherical harmonic gravity fields.
C++ Library with command-line tools and Python bindings.

![Earth Gravity Anomalies](https://github.com/DavidMSCode/Spherical-Harmonic-Gravity/blob/main/examples/Anomalies1800.png)

## Features

- Calculate gravitational acceleration and potential from spherical harmonic coefficients
- Support for arbitrary gravity models (EGM2008, custom fields, planetary bodies)
- Maximum degree ~1900 (unmodified forward column method limitation)
- C++17 library with Python bindings
- Demo programs and visualization examples included

## Quick Start

### C++ Build
```bash
# Clone repository
git clone https://github.com/DavidMSCode/Spherical-Harmonic-Gravity.git
cd Spherical-Harmonic-Gravity

# Build with CMake
mkdir build && cd build
cmake ..
make

# Run C++ examples
./GravityCalc      # Gravity field calculator
./PotentialCalc    # Gravitational potential calculator
./Demo             # Demo program that calculates gravity and potential at sample locations
```

### Python Build
```bash
# Clone repository
git clone https://github.com/DavidMSCode/Spherical-Harmonic-Gravity.git
cd Spherical-Harmonic-Gravity

# Install Python dependencies
pip install numpy matplotlib cartopy pybind11

# Build and install Python module
python setup.py build
python setup.py install

# Or use pip for development install
pip install -e .
```

### EGM2008 Coefficients Setup
```bash
# Download EGM2008 coefficients from NGA (104MB)
wget https://earth-info.nga.mil/php/download.php?file=egm-08spherical
# Or manually visit https://earth-info.nga.mil and download the EGM2008 spherical harmonics file.

# Extract and place coefficients file EGM2008_to2190_TideFree in the project directory you plan to run from.
# File format: ASCII text with degree, order, C_nm, S_nm coefficients
```

### Geopotential artifact bins
The preferred C++ API is model-based and loads coefficients from the shared artifact bins.
Supported model IDs in the shared artifact release:
- `EGM2008`
- `EGM96`
- `GRGM1200A`
- `GMM3`

The default functions for potential and acceleration take a model object as the first argument, and the library provides helpers to get models by ID or body name. For example:
- `SHG::get_model("EGM2008")`
- `SHG::default_model_for_body("Earth")`
- `SHG::potential("EGM2008", r, phi, lambda)`
- `SHG::acceleration("EGM2008", r, phi, lambda)`
- `SHG::potential(model, r, phi, lambda)`
- `SHG::acceleration(model, r, phi, lambda)`

Example:
```cpp
#include <array>

SHG::SHGModel model = SHG::get_model(SHG::default_model_for_body("Earth"));
double U = SHG::potential(model, r, phi, lambda);
std::array<double, 3> g = SHG::acceleration(model, r, phi, lambda);
```

You can set the artifact location with `SHG::set_model_bins_directory(...)`, or
by setting `GEOPOTENTIAL_MODEL_BINS_DIR` before running your program.

By default, CMake fetches the bins during the build so the library and tests can
use them immediately. If you need an offline workflow, point CMake at a local
tarball with `-DSHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE=/path/to/geopotential_bins.tar.gz`
and it will unpack that instead of trying to download them.

If you want to suppress automatic fetching entirely, disable it with
`-DSHG_DOWNLOAD_GEOPOTENTIAL_ARTIFACTS=OFF` and build the
`download_geopotential_artifacts` target only when you want to populate the cache.

### Compatibility / legacy usage
```python
# Compute gravitational acceleration at a point
import pyshg
r, phi, lam = 6371000, 0.5, 1.0  # radius (m), geocentric latitude (rad), longitude (rad)
g_vector = pyshg.g_EGM2008(r, phi, lam, max_degree=360)

# The direct coefficient-path workflow remains available for existing scripts.

# Run demonstration examples
python examples/earth_gravity_visualization.py
python examples/load_gravity_data_demo.py
```

## Examples

```bash
python examples/earth_gravity_visualization.py
python examples/load_gravity_data_demo.py
```

## Technical Notes

- Uses unmodified forward column method for computing the associate Legendre functions (order limit ~1900)
- Supports spherical coordinates (r, theta, phi) in geocentric frame and cartesian (x, y, z) coordinates
- Compatible with EGM2008 and custom gravity field coefficients

## References

- Holmes, S., Featherstone, W. A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions. Journal of Geodesy 76, 279–299 (2002). https://doi.org/10.1007/s00190-002-0216-2
- "Fundamentals of Astrodynamics and Applications" by David A. Vallado, 4th edition, Microcosm Press, 2013. 

## License

MIT License - See LICENSE file for details.
