.PHONY: all clean cython python

# Install the python packages needed to run SOAPFAST

PIPOPTS=--user

all: python cython

python:
	pip install $(PIPOPTS) numpy==1.16.2 ase==3.16.2
	pip install $(PIPOPTS) scipy cython sympy

cython:
	cd utils;python setup_cython.py build_ext --inplace

clean:
	cd utils;rm -rf *.so build
