#!/bin/bash

FC=gfortran
FCFLAGS="-fbacktrace -Wall -Wtabs"
FCDEBUG="-Og -g ${FCFLAGS} -Wextra -fcheck=all"
LPFLAGS="-llapack"

PROGRAM=cmplx_lyap

[[ -f "${PROGRAM}" ]] && rm *.mod *.o ${PROGRAM}
[[ "$1" -eq "d" ]] && FCF=${FCDEBUG} || FCF=${FCFLAGS}

gfortran -c mod_prec.f90 ${FCF}
gfortran -c mod_proc.f90 ${FCF}
gfortran -c mod_lpck.f90 ${FCF} ${LPFLAGS}
gfortran -o ${PROGRAM} prog_lyap.f90 mod_prec.o mod_proc.o mod_lpck.o ${LPFLAGS}
