#ifndef UTILS_H
#define UTILS_H

#define NO_IMPORT_ARRAY
#include "common.h"

int not_float64_vector(PyArrayObject *);
npy_intp ndarray_dim(PyArrayObject *);
double get_double_at(PyArrayObject *, npy_intp);
double* reverse_vector(const double*, npy_intp);

#endif