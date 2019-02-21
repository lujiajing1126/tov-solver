#include "utils.h"

/* ==== Check that PyArrayObject is a double (Float) type and a vector ==============
    return 1 if an error and raise exception */ 
int not_float64_vector(PyArrayObject *vec)  {
	if (PyArray_DESCR(vec)->type_num != NPY_FLOAT64 || PyArray_NDIM(vec) != 1)  {
		PyErr_SetString(PyExc_ValueError,
			"In not_doublevector: array must be of type Float and 1 dimensional (n).");
		return 1;  }
	return 0;
}

npy_intp ndarray_dim(PyArrayObject *vec) {
    return PyArray_SIZE(vec);
}

double get_double_at(PyArrayObject *array, npy_intp idx) {
	PyObject* val = PyArray_GETITEM(array, PyArray_GETPTR1(array, idx));
	if (val == NULL) {
		return 0.0;
	}
	if (PyFloat_Check(val)) {
		return PyFloat_AsDouble(val);
	} else {
		return 0.0;
	}
}

double *reverse_vector(const double *pointer, npy_intp n) {
	npy_intp c = n - 1;
	double *s = (double*)PyMem_Malloc(sizeof(double)*n);
	if ( s == NULL ) {
		exit(EXIT_FAILURE);
	}
	for (;c>=0;c--) {
		*(s+n-1-c) = *(pointer+c);
	}
	return s;
}