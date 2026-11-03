#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP drw_estimate(SEXP,SEXP,SEXP,SEXP);
extern SEXP drw_outcome_likelihood(SEXP,SEXP,SEXP);
extern SEXP drw_cox_likelihood(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"drw_estimate", (DL_FUNC) &drw_estimate, 4},
    {"drw_outcome_likelihood", (DL_FUNC) &drw_outcome_likelihood, 3},
    {"drw_cox_likelihood", (DL_FUNC) &drw_cox_likelihood, 9},
    {NULL, NULL, 0}
};

void attribute_visible R_init_DRwinStats(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
