#include "solver.h"

// methods declare
static PyObject* tov_solve(PyObject*, PyObject*);
static Pair* evaluate(double, double, double, double, gsl_spline *, gsl_interp_accel*, gsl_spline *, gsl_interp_accel*);

// method table
static PyMethodDef TOVMethods[] = {
    {"solve",  (PyCFunction)tov_solve, METH_VARARGS,
     "Entry for solving the equation."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

// module definition
static struct PyModuleDef tov_module = {
    PyModuleDef_HEAD_INIT,
    "tov",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    TOVMethods
};

// Init Module
PyMODINIT_FUNC PyInit_tov(void)  {
    PyObject *m;
    m = PyModule_Create(&tov_module);

    if (m == NULL) {
        return NULL;
    }
    
    _import_array();

    return m;
}

/**
 * Entry for solving TOV Equations.
 * 
 * Python signature:
 * def solve(rho: np.ndarray, press: np.ndarray, eps: np.ndarray)
 */
static PyObject* tov_solve(PyObject *dummy, PyObject *args) {
    PyArrayObject *rho, *press, *eps;
    gsl_spline *spline_e_p;
    gsl_spline *spline_rho_p;
    gsl_interp_accel *accel_e_p;
    gsl_interp_accel *accel_rho_p;

    if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &rho, &PyArray_Type, &press, &PyArray_Type, &eps))
        return NULL;

    if (NULL == rho) return NULL;
	if (NULL == press) return NULL;
    if (NULL == eps) return NULL;
    
    // ensure the dtype of the ndarray(s) is float64
    if (not_float64_vector(rho)) return NULL;
    if (not_float64_vector(press)) return NULL;
    if (not_float64_vector(eps)) return NULL;

    // ensure the three vector have the same size
    if (ndarray_dim(rho) != ndarray_dim(press)) return NULL;
    if (ndarray_dim(rho) != ndarray_dim(eps)) return NULL;

    npy_intp size = ndarray_dim(rho);

    // Pass in NPY_MAXDIMS for axis in order to achieve the same effect 
    // that is obtained by passing in axis = None in Python (treating the array as a 1-d array).
    // https://docs.scipy.org/doc/numpy-1.16.1/reference/c-api.array.html#calculation
    PyObject *min_pressure_obj = NULL;
    double min_pressure = 0.0;
    min_pressure_obj = PyArray_Min(press, NPY_MAXDIMS, NULL);
    if (min_pressure_obj == NULL) {
        printf("cannot get min value of pressure");
        return NULL;
    }
    // check return object is scalar
    if (PyArray_IsScalar(min_pressure_obj, Double)) {
        PyArray_ScalarAsCtype(min_pressure_obj, &min_pressure);
        printf("Minimum Pressure: %.2f\n", min_pressure);
        const double *eps_c = reverse_vector((double*) PyArray_DATA(eps), size);
        const double *press_c = reverse_vector((double*) PyArray_DATA(press), size);
        const double *rho_c = reverse_vector((double*) PyArray_DATA(rho), size);
        spline_e_p = gsl_spline_alloc(gsl_interp_akima, size);
        accel_e_p = gsl_interp_accel_alloc();
        spline_rho_p = gsl_spline_alloc(gsl_interp_akima, size);
        accel_rho_p = gsl_interp_accel_alloc();
        gsl_spline_init(spline_e_p, press_c, eps_c, size);
        gsl_spline_init(spline_rho_p, press_c, rho_c, size);
        for (npy_intp i = 0; i < size; i++) {
            double rho0 = get_double_at(rho, i);
            if (rho0 < 0.10) {
                continue;
            }
            double P0 = get_double_at(press, i), MG0 = 0.0, MB0 = 0.0, R = 1e-10, DR = 1e-5;
            for (int j = 0; j < 1000000; j++) {
                Pair *retValue = evaluate(R, P0, MG0, min_pressure, spline_e_p, accel_e_p, spline_rho_p, accel_rho_p);
                double DP1 = retValue->a, DM1 = retValue->b, DMB1 = retValue->c;
                PyMem_Free(retValue);

                if (j > 1) {
                    DR = 0.01/(DM1/MG0 - DP1/P0);
                }

                R = R + DR;

                double DR2 = DR/2.0;

                // -> K2
                Pair *retValue2 = evaluate(R+DR2, P0+DR2*DP1, MG0+DR2*DM1, min_pressure, spline_e_p, accel_e_p, spline_rho_p, accel_rho_p);
                double DP2 = retValue2->a, DM2 = retValue2->b, DMB2 = retValue2->c; 
                PyMem_Free(retValue2);
                // -> K3
                Pair *retValue3 = evaluate(R+DR2, P0+DR2*DP2, MG0+DR2*DM2, min_pressure, spline_e_p, accel_e_p, spline_rho_p, accel_rho_p);
                double DP3 = retValue3->a, DM3 = retValue3->b, DMB3 = retValue3->c; 
                PyMem_Free(retValue3);
                // -> K4
                Pair *retValue4 = evaluate(R+DR,  P0+DR*DP3,  MG0+DR*DM3, min_pressure, spline_e_p, accel_e_p, spline_rho_p, accel_rho_p);
                double DP4 = retValue4->a, DM4 = retValue4->b, DMB4 = retValue4->c;
                PyMem_Free(retValue4);

                //...NEW P AND M...
                P0  = P0  + DR * (DP1  + 2*DP2  + 2*DP3  + DP4 ) / 6;
                MG0 = MG0 + DR * (DM1  + 2*DM2  + 2*DM3  + DM4 ) / 6;
                MB0 = MB0 + DR * (DMB1 + 2*DMB2 + 2*DMB3 + DMB4) / 6;

                if (P0 < min_pressure) {
                    break;
                }
            }
            printf("%ld %.4f %.4f %.4f %.4f %.4f %.4f\n", i, MG0/MS, R/1e5, get_double_at(eps,i)/CE, get_double_at(press,i)/CP, MB0/MS, rho0);
        }
        gsl_spline_free(spline_e_p);
        gsl_spline_free(spline_rho_p);
    } else {
        return NULL;
    }

    Py_RETURN_NONE;
}

static Pair* evaluate(double R0, double P0, double MG0, double minimum_p, gsl_spline *spline_e_p, gsl_interp_accel *accel_e_p, gsl_spline *spline_rho_p, gsl_interp_accel *accel_rho_p) {    
    // Use PyMem_RawMalloc
    // https://docs.python.org/3/c-api/memory.html#c.PyMem_RawMalloc
    Pair* pair = (Pair*) PyMem_RawMalloc(sizeof(Pair));
    if (P0 < minimum_p) {
        pair->a = 0.0;
        pair->b = 0.0;
        pair->c = 0.0;
        return pair;
    }
    // 10.1103/PhysRevC.100.054335 Eq. (17), (18), (19)
    double E = gsl_spline_eval(spline_e_p, P0, accel_e_p);
    double D1 = - GC * (E + P0/C2) * (MG0/R0 + PI4*R0*R0*P0/C2) / (R0 - 2*MG0*GC/C2);
    double D2 = PI4 * R0*R0 * E;
    double RHO = gsl_spline_eval(spline_rho_p, P0, accel_rho_p);
    double D3 = PI4 * R0*R0 * RHO * MN / sqrt(1.0 - 2*GC*MG0/R0/C2);
    pair->a = D1;
    pair->b = D2;
    pair->c = D3;
    return pair;
}