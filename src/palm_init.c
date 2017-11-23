#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _palm_buffer_keep(SEXP, SEXP, SEXP);
extern SEXP _palm_pbc_distances(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_palm_buffer_keep",   (DL_FUNC) &_palm_buffer_keep,   3},
    {"_palm_pbc_distances", (DL_FUNC) &_palm_pbc_distances, 2},
    {NULL, NULL, 0}
};

void R_init_palm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
