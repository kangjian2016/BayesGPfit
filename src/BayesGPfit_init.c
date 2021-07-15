#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .C calls */
extern void Wrapper_R_GP_eigen_funcs(void *, void *, void *, void *, void *, void *, void *, void *);
extern void Wrapper_R_GP_eigen_funcs_orth(void *, void *, void *, void *, void *, void *, void *, void *);

//extern void Wrapper_R_GP_eigen_funcs(double* , int* double*, int*, int*, int*, double*, double*);
//extern void Wrapper_R_GP_eigen_funcs_orth(double* , int* double*, int*, int*, int*, double*, double*);


static const R_CMethodDef CEntries[] = {
  {"Wrapper_R_GP_eigen_funcs",      (DL_FUNC) &Wrapper_R_GP_eigen_funcs,      8},
  {"Wrapper_R_GP_eigen_funcs_orth", (DL_FUNC) &Wrapper_R_GP_eigen_funcs_orth, 8},
  {NULL, NULL, 0}
};

void R_init_BayesGPfit(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
