from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

ext_modules = [
        Extension(
        "fourier_integrals",
        ["fourier_integrals.pyx"],
        libraries=["m"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    ),
        Extension(
        "realspace_potential",
        ["realspace_potential.pyx"],
        libraries=["m"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    ),
        Extension(
        "fourier_integrals_MT2D",
        ["fourier_integrals_MT2D.pyx"],
        libraries=["m"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
]

setup(
    name='Cython routines for LODE power spectrum evaluation',
    ext_modules=cythonize(ext_modules),
)
