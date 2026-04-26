#ifndef GLMBAYES_GETREGISTEREDNAMESPACE_H
#define GLMBAYES_GETREGISTEREDNAMESPACE_H

#include <Rversion.h>

/*
 * Package-local compatibility shim for Rcpp Function.h namespace lookup.
 * Only active on R toolchains where Rcpp may select R_getRegisteredNamespace
 * but headers / API snapshot are not aligned.
 */
#if R_VERSION < R_Version(4,6,0) || R_SVN_REVISION < 89746
#define GLMBAYES_NEED_GETREGISTEREDNAMESPACE_SHIM 1
struct SEXPREC;
typedef SEXPREC* SEXP;
extern "C" SEXP glmbayes_getRegisteredNamespace(const char *name);
#define R_getRegisteredNamespace glmbayes_getRegisteredNamespace
#endif

#endif
