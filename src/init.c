#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _MagmaClustR_cpp_dist(SEXP, SEXP);
extern SEXP _MagmaClustR_cpp_noise(SEXP, SEXP, SEXP);
extern SEXP _MagmaClustR_cpp_perio(SEXP, SEXP, SEXP);
extern SEXP _MagmaClustR_cpp_perio_deriv(SEXP, SEXP, SEXP);
extern SEXP _MagmaClustR_cpp_prod(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_MagmaClustR_cpp_dist",        (DL_FUNC) &_MagmaClustR_cpp_dist,        2},
    {"_MagmaClustR_cpp_noise",       (DL_FUNC) &_MagmaClustR_cpp_noise,       3},
    {"_MagmaClustR_cpp_perio",       (DL_FUNC) &_MagmaClustR_cpp_perio,       3},
    {"_MagmaClustR_cpp_perio_deriv", (DL_FUNC) &_MagmaClustR_cpp_perio_deriv, 3},
    {"_MagmaClustR_cpp_prod",        (DL_FUNC) &_MagmaClustR_cpp_prod,        2},
    {NULL, NULL, 0}
};

void R_init_MagmaClustR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
