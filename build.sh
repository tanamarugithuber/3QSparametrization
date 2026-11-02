#!/bin/bash
# FRLDM + Shell + pairing (main.exe). Needs FFTW3 and LAPACK.
# FFTW3 の include 先が標準以外なら FFTW_INC (と FFTW_LIB) で指定:
#   FFTW_INC=/opt/fftw/include FFTW_LIB=/opt/fftw/lib ./build.sh
set -e

if [ -z "$FFTW_INC" ]; then
    FFTW_INC=$(pkg-config --variable=includedir fftw3 2>/dev/null || true)
    [ -n "$FFTW_INC" ] || FFTW_INC=/usr/include
fi
if [ ! -f "$FFTW_INC/fftw3.f03" ]; then
    echo "error: $FFTW_INC/fftw3.f03 が見つかりません (FFTW3 の開発用パッケージが必要)" >&2
    echo "  Ubuntu/Debian: sudo apt install libfftw3-dev liblapack-dev" >&2
    echo "  Fedora/RHEL  : sudo dnf install fftw-devel lapack-devel" >&2
    echo "  macOS        : brew install fftw lapack" >&2
    echo "  MSYS2        : pacman -S mingw-w64-x86_64-fftw mingw-w64-x86_64-lapack" >&2
    echo "  別の場所にある場合: FFTW_INC=<fftw3.f03 のあるディレクトリ> ./build.sh" >&2
    exit 1
fi
LIBS="-lfftw3 -llapack"
[ -n "$FFTW_LIB" ] && LIBS="-L$FFTW_LIB $LIBS"

rm -f *.mod
rm -f *.exe
rm -f *.o

export OMP_NUM_THREADS=4
FC="gfortran -O2 -fopenmp"

$FC -c constant_mod.f90
$FC -c nucleus_mod.f90
$FC -c shape_mod.f90
$FC -c qs3_shape_mod.f90
$FC -c grid_mod.f90
$FC -I"$FFTW_INC" -c calculate_potential_hockney.f90
$FC -c frldm_mod.f90
$FC -c micro_constant_mod.f90
$FC -c ho_basis_mod.f90
$FC -c potential_gh_mod.f90
$FC -c shell_bcs_mod.f90
$FC -c zero_point_mod.f90
$FC -c main.f90

$FC \
    constant_mod.o \
    nucleus_mod.o \
    shape_mod.o \
    qs3_shape_mod.o \
    grid_mod.o \
    calculate_potential_hockney.o \
    frldm_mod.o \
    micro_constant_mod.o \
    ho_basis_mod.o \
    potential_gh_mod.o \
    shell_bcs_mod.o \
    zero_point_mod.o \
    main.o \
    $LIBS \
    -o main.exe
