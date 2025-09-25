#!/usr/bin/env python3
"""
Setup script for Python bindings of Spherical Harmonic Gravity (SHG) library
"""

try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext, ParallelCompile
    from pybind11.setup_helpers import Pybind11Extension
    from setuptools import setup, Extension
    from pybind11 import get_cmake_dir
    import pybind11
except ImportError:
    import subprocess
    import sys
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pybind11'])
    from pybind11.setup_helpers import Pybind11Extension, build_ext, ParallelCompile
    from pybind11.setup_helpers import Pybind11Extension
    from setuptools import setup, Extension
    from pybind11 import get_cmake_dir
    import pybind11

__version__ = "1.0.0"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=14/17/20, or even cxx_std="c++17"
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

ext_modules = [
    Pybind11Extension(
        "pyshg",
        [
            "python/shg_python.cpp",
            "src/SHG.cpp",
        ],
        include_dirs=[
            # Path to pybind11 headers
            pybind11.get_include(),
            "include",
        ],
        cxx_std=17,
        define_macros=[("VERSION_INFO", '"{}"'.format(__version__))],
    ),
]

setup(
    name="pyshg",
    version=__version__,
    author="DavidMSCode", 
    author_email="",
    url="https://github.com/DavidMSCode/Spherical-Harmonic-Gravity",
    description="Python bindings for Spherical Harmonic Gravity (SHG) library",
    long_description="""
    Python bindings for a C++ library that computes gravitational potential and 
    acceleration using normalized spherical harmonic coefficients. The library 
    supports high-degree spherical harmonic gravity models such as EGM2008.
    
    Features:
    - Compute gravitational acceleration and potential
    - Support for both geocentric and Cartesian coordinates
    - Efficient normalized Associated Legendre Function computation
    - Binary and text coefficient file I/O
    - Template-based design for flexibility
    """,
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.6",
    setup_requires=["pybind11>=2.5.0"],
    install_requires=[
        "numpy>=1.15.0",
        "pybind11>=2.5.0",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research", 
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
)