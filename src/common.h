#ifndef COMMON_H
#define COMMON_H
// include python core headers and module-scope headers
#include <Python.h>
// Not use deprecated API in NumPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// Solve the initialization problem,
// Use unique py_array symbol
#define PY_ARRAY_UNIQUE_SYMBOL tov_ARRAY_API
#include <numpy/arrayobject.h>

#include "const.h"

#endif