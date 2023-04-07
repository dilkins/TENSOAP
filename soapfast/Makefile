.PHONY: all clean cython python LODE

F2PY=f2py3
PIP=pip3
PYTH=python3
CC=gcc

# Install the python packages needed to run SOAPFAST

PIPOPTS=--user

all: cython

LODE: gvectors.so phasecomb.so gcontra.so lebedev.so gausslegendre.so nearfield.so nearfield_ewald.so neighlist_ewald.so potrefewald.so potrefewaldmt2d.so lode-cython

gvectors.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' gvectors.f90 -m gvectors;mv gvectors.*.so gvectors.so

phasecomb.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' phasecomb.f90 -m phasecomb --f90flags="-fopenmp" -lgomp;mv phasecomb.*.so phasecomb.so

nearfield.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' nearfield.f90 -m nearfield --f90flags="-fopenmp" -lgomp;mv nearfield.*.so nearfield.so

nearfield_ewald.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' nearfield_ewald.f90 -m nearfield_ewald --f90flags="-fopenmp" -lgomp;mv nearfield_ewald.*.so nearfield_ewald.so

neighlist_ewald.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' neighlist_ewald.f90 -m neighlist_ewald --f90flags="-fopenmp" -lgomp;mv neighlist_ewald.*.so neighlist_ewald.so

gcontra.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' gcontra.f90 -m gcontra --f90flags="-fopenmp" -lgomp;mv gcontra.*.so gcontra.so

lebedev.so:
	cd utils/LODE; $(CC) -fPIC -shared -o lebedev.so Lebedev-Laikov.c

gausslegendre.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' gausslegendre.f90 -m gausslegendre;mv gausslegendre.*.so gausslegendre.so

potrefewald.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' potrefewald.f90 -m potrefewald;mv potrefewald.*.so potrefewald.so

potrefewaldmt2d.so:
	cd utils/LODE;$(F2PY) -c --opt='-O2' potrefewaldmt2d.f90 -m potrefewaldmt2d;mv potrefewaldmt2d.*.so potrefewaldmt2d.so

lode-cython:
	cd utils/LODE;$(PYTH) setup_cython.py build_ext --inplace
	@cd utils/LODE;mv fourier_integrals.*.so fourier_integrals.so
	@cd utils/LODE;mv fourier_integrals_MT2D.*.so fourier_integrals_MT2D.so
	@cd utils/LODE;mv realspace_potential.*.so realspace_potential.so

install-deps:
	$(PIP) install $(PIPOPTS) numpy ase scipy cython sympy

cython: #install-deps
	cd utils;$(PYTH) setup_cython.py build_ext --inplace
	@cd utils;mv initsoap.*.so initsoap.so;mv build_kernel.*.so build_kernel.so

clean:
	cd utils;rm -rf *.so build;cd LODE;rm -rf *.so build
