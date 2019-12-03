#ifndef SOLVER_H
#define SOLVER_H

#if PY_MAJOR_VERSION >= 3
#define PY3K
#endif

// include user-defined headers
#include "common.h"
#include "utils.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

typedef struct _TOV_Tuple {
    double a;
    double b;
    double c;
} Pair;

#endif