#!/bin/bash
# 3QS fission saddle points: pes3qs.exe (potential-energy surface) and saddle3qs.exe (immersion)
# Needs FFTW3 and LAPACK (as build.sh). FFTW_INC / FFTW_LIB as in build.sh.
set -e

if [ -z "$FFTW_INC" ]; then
    FFTW_INC=$(pkg-config --variable=includedir fftw3 2>/dev/null || true)
    [ -n "$FFTW_INC" ] || FFTW_INC=/usr/include
fi
if [ ! -f "$FFTW_INC/fftw3.f03" ]; then
    echo "error: $FFTW_INC/fftw3.f03 が見つかりません (FFTW3 の開発用パッケージが必要, build.sh を参照)" >&2
    exit 1
fi
LIBS="-lfftw3 -llapack"
[ -n "$FFTW_LIB" ] && LIBS="-L$FFTW_LIB $LIBS"

FC="gfortran -O2 -fopenmp"
B=build_3qs
mkdir -p $B
SRC="constant_mod nucleus_mod shape_mod qs3_shape_mod qs3_coord_mod grid_mod calculate_potential_hockney
     frldm_mod micro_constant_mod ho_basis_mod potential_gh_mod shell_bcs_mod immersion_mod"
OBJ=""
for f in $SRC; do
    $FC -I"$FFTW_INC" -J$B -c $f.f90 -o $B/$f.o
    OBJ="$OBJ $B/$f.o"
done

$FC -J$B pes3qs.f90 $OBJ $LIBS -o pes3qs.exe
$FC -J$B saddle3qs.f90 $B/immersion_mod.o -o saddle3qs.exe

# tests
$FC -J$B test_qs3.f90 $B/constant_mod.o $B/shape_mod.o $B/qs3_shape_mod.o $B/qs3_coord_mod.o -o test_qs3.exe
$FC -J$B test_frldm_3qs.f90 $OBJ $LIBS -o test_frldm_3qs.exe
$FC -J$B test_immersion.f90 $B/immersion_mod.o -o test_immersion.exe
echo "built: pes3qs.exe saddle3qs.exe test_qs3.exe test_frldm_3qs.exe test_immersion.exe"
