#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(mybnorm)(void *, void *, void *, void *, void *);
extern void F77_NAME(mytnorm)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(myunorm)(void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"mybnorm", (DL_FUNC) &F77_NAME(mybnorm), 5},
    {"mytnorm", (DL_FUNC) &F77_NAME(mytnorm), 6},
    {"myunorm", (DL_FUNC) &F77_NAME(myunorm), 3},
    {NULL, NULL, 0}
};

void R_init_mhurdle(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
