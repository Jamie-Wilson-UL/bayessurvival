#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(lddpsurvival)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hpd)(void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"lddpsurvival", (DL_FUNC) &F77_NAME(lddpsurvival), 64},
    {"hpd",          (DL_FUNC) &F77_NAME(hpd),           5},
    {NULL, NULL, 0}
};

#ifdef _WIN32
# include <fcntl.h>
#endif

void R_init_bayessurvival(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);

#ifdef _WIN32
    /* gfortran initialization sets these to _O_BINARY */
    setmode(1, _O_TEXT); /* stdout */
    setmode(2, _O_TEXT); /* stderr */
#endif
}
