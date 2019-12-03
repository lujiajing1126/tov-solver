#ifndef CONST_H
#define CONST_H

#include <numpy/npy_math.h>

#define PI4 (NPY_PI * 4.0)
// Gravity Constant in CGS
// CM^3/GM S^2
#define GC 6.6742e-8
#define VC 2.99792458e10
#define CE 1.78266e12
#define MS 1.9891e33
#define C2 (VC * VC)
#define CP (CE * C2)
// Nucleon Mass in CGS
#define MN 1.67e15

#endif