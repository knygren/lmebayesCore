# Synchronize installed Rcpp Function.h with the local R headers.
# If Function.h appears "frozen" to one path that conflicts with this R's
# R_VERSION / R_SVN_REVISION branch, replace it with the upstream conditional
# header so compile-time preprocessor selection chooses the right path.

read_r_svn_revision_h <- function() {
  rh <- file.path(R.home("include"), "Rversion.h")
  if (!file.exists(rh)) {
    return(NA_integer_)
  }
  z <- readLines(rh, warn = FALSE, encoding = "UTF-8")
  m <- grep("^#define[[:space:]]+R_SVN_REVISION[[:space:]]+", z, value = TRUE)
  if (!length(m)) {
    return(NA_integer_)
  }
  rest <- sub("^#define[[:space:]]+R_SVN_REVISION[[:space:]]+", "", m[1L])
  rest <- sub("/\\*.*$", "", rest)
  rest <- gsub("^[[:space:]]+|[[:space:]]+$", "", rest)
  suppressWarnings(as.integer(rest))
}

expected_branch <- function() {
  rv <- getRversion()
  svn_h <- read_r_svn_revision_h()
  if (rv < "4.5.0") {
    return(1L)
  }
  if (rv < "4.6.0" || (!is.na(svn_h) && svn_h < 89746L)) {
    return(2L)
  }
  3L
}

is_universal_header <- function(txt) {
  grepl("#if[[:space:]]+R_VERSION[[:space:]]*<[[:space:]]*R_Version\\(4,5,0\\)", txt, perl = TRUE) &&
    grepl("#elif[[:space:]]+R_VERSION[[:space:]]*<[[:space:]]*R_Version\\(4,6,0\\)[[:space:]]*\\|\\|[[:space:]]*R_SVN_REVISION[[:space:]]*<[[:space:]]*89746", txt, perl = TRUE) &&
    grepl("R_getRegisteredNamespace\\(", txt, perl = TRUE)
}

detected_frozen_branch <- function(txt) {
  has_b1 <- grepl("Rf_findVarInFrame\\([^)]*R_NamespaceRegistry", txt, perl = TRUE)
  has_b2 <- grepl("R_getVarEx\\([^)]*R_NamespaceRegistry", txt, perl = TRUE)
  has_b3 <- grepl("R_getRegisteredNamespace\\(", txt, perl = TRUE)
  n <- sum(c(has_b1, has_b2, has_b3))
  if (n != 1L) {
    return(NA_integer_)
  }
  if (has_b1) return(1L)
  if (has_b2) return(2L)
  3L
}

full_function_h <- c(
  "// Function.h: Rcpp R/C++ interface class library -- functions (also primitives and builtins)",
  "//",
  "// Copyright (C) 2010 - 2025  Dirk Eddelbuettel and Romain Francois",
  "// Copyright (C) 2026         Dirk Eddelbuettel, Romain Francois and Inaki Ucar",
  "//",
  "// This file is part of Rcpp.",
  "//",
  "// Rcpp is free software: you can redistribute it and/or modify it",
  "// under the terms of the GNU General Public License as published by",
  "// the Free Software Foundation, either version 2 of the License, or",
  "// (at your option) any later version.",
  "//",
  "// Rcpp is distributed in the hope that it will be useful, but",
  "// WITHOUT ANY WARRANTY; without even the implied warranty of",
  "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the",
  "// GNU General Public License for more details.",
  "//",
  "// You should have received a copy of the GNU General Public License",
  "// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.",
  "",
  "#ifndef Rcpp_Function_h",
  "#define Rcpp_Function_h",
  "",
  "#include <RcppCommon.h>",
  "",
  "#include <Rcpp/grow.h>",
  "",
  "namespace Rcpp{",
  "",
  "    /**",
  "     * functions",
  "     */",
  "    RCPP_API_CLASS(Function_Impl) {",
  "    public:",
  "",
  "        RCPP_GENERATE_CTOR_ASSIGN(Function_Impl)",
  "",
  "        Function_Impl(SEXP x){",
  "            switch( TYPEOF(x) ){",
  "            case CLOSXP:",
  "            case SPECIALSXP:",
  "            case BUILTINSXP:",
  "                Storage::set__(x);",
  "                break;",
  "            default:                        // #nocov start",
  "                const char* fmt = \"Cannot convert object to a function: \"",
  "                                  \"[type=%s; target=CLOSXP, SPECIALSXP, or \"",
  "                                  \"BUILTINSXP].\";",
  "                throw not_compatible(fmt, Rf_type2char(TYPEOF(x)));",
  "            }                               // #nocov end",
  "        }",
  "",
  "        /**",
  "         * Finds a function. By default, searches from the global environment",
  "         *",
  "         * @param name name of the function",
  "         * @param env an environment where to search the function",
  "         * @param ns name of the namespace in which to search the function",
  "         */",
  "        Function_Impl(const std::string& name) {",
  "            get_function(name, R_GlobalEnv);",
  "        }",
  "",
  "        Function_Impl(const std::string& name, const SEXP env) {",
  "            if (!Rf_isEnvironment(env)) {",
  "                stop(\"env is not an environment\");",
  "            }",
  "            get_function(name, env);",
  "        }",
  "",
  "        Function_Impl(const std::string& name, const std::string& ns) {",
  "#if R_VERSION < R_Version(4,5,0)",
  "            // before R 4.5.0 we would use Rf_findVarInFrame",
  "            Shield<SEXP> env(Rf_findVarInFrame(R_NamespaceRegistry, Rf_install(ns.c_str())));",
  "            if (env == R_UnboundValue)",
  "                stop(\"there is no namespace called \\\"%s\\\"\", ns);",
  "#elif R_VERSION < R_Version(4,6,0) || R_SVN_REVISION < 89746",
  "            // during R 4.5.* and before final R 4.6.0 we could use R_getVarEx",
  "            // along with R_NamespaceRegistry but avoid R_UnboundValue",
  "            Shield<SEXP> env(R_getVarEx(Rf_install(ns.c_str()), R_NamespaceRegistry, FALSE, R_NilValue));",
  "            if (env == R_NilValue)",
  "                stop(\"there is no namespace called \\\"%s\\\"\", ns);",
  "#else",
  "            // late R 4.6.0 development got us R_getRegisteredNamespace",
  "            Shield<SEXP> env(R_getRegisteredNamespace(ns.c_str()));",
  "            if (env == R_NilValue)",
  "                stop(\"there is no namespace called \\\"%s\\\"\", ns);",
  "#endif",
  "            get_function(name, env);",
  "        }",
  "",
  "        SEXP operator()() const {",
  "            Shield<SEXP> call(Rf_lang1(Storage::get__()));",
  "            return Rcpp_fast_eval(call, R_GlobalEnv);",
  "        }",
  "",
  "        template <typename... T>",
  "        SEXP operator()(const T&... args) const {",
  "            return invoke(pairlist(args...), R_GlobalEnv);",
  "        }",
  "",
  "        /**",
  "         * Returns the environment of this function",
  "         */",
  "        SEXP environment() const {",
  "            SEXP fun = Storage::get__() ;",
  "            if( TYPEOF(fun) != CLOSXP ) {",
  "                throw not_a_closure(Rf_type2char(TYPEOF(fun)));",
  "            }",
  "            #if (defined(R_VERSION) && R_VERSION >= R_Version(4,5,0))",
  "            return R_ClosureEnv(fun);",
  "            #else",
  "            return CLOENV(fun);",
  "            #endif",
  "        }",
  "",
  "        /**",
  "         * Returns the body of the function",
  "         */",
  "        SEXP body() const {",
  "            return BODY( Storage::get__() ) ;",
  "        }",
  "",
  "        void update(SEXP){}",
  "",
  "",
  "    private:",
  "        void get_function(const std::string& name, const SEXP env) {",
  "            SEXP nameSym = Rf_install( name.c_str() );    // cannot be gc()'ed  once in symbol table",
  "            Shield<SEXP> x( Rf_findFun( nameSym, env ) ) ;",
  "            Storage::set__(x) ;",
  "        }",
  "",
  "        SEXP invoke(SEXP args_, SEXP env) const {",
  "            Shield<SEXP> args(args_);",
  "            Shield<SEXP> call(Rcpp_lcons(Storage::get__(), args));",
  "            SEXP out = Rcpp_fast_eval(call, env);",
  "            return out;",
  "        }",
  "",
  "    };",
  "",
  "    typedef Function_Impl<PreserveStorage> Function ;",
  "",
  "} // namespace Rcpp",
  "",
  "#endif"
)

lib <- Sys.getenv("GLMBAYES_RCPP_LIB", "")
if (!nzchar(lib)) lib <- Sys.getenv("R_LIBS_USER", "")
if (!nzchar(lib)) lib <- .libPaths()[1L]
fh <- file.path(lib, "Rcpp", "include", "Rcpp", "Function.h")

if (!file.exists(fh)) {
  message("patch_rcpp_function_h: missing ", fh, " - skip")
  quit(status = 0L)
}

lines <- readLines(fh, warn = FALSE, encoding = "UTF-8")
txt <- paste(lines, collapse = "\n")
exp_branch <- expected_branch()

if (is_universal_header(txt)) {
  message("patch_rcpp_function_h: Function.h already has full conditional branch logic - skip")
  quit(status = 0L)
}

det_branch <- detected_frozen_branch(txt)
if (!is.na(det_branch) && det_branch == exp_branch) {
  message(
    "patch_rcpp_function_h: frozen Function.h appears consistent with this R (detected branch ",
    det_branch, ", expected ", exp_branch, ") - skip"
  )
  quit(status = 0L)
}

backup <- paste0(fh, ".glmbayes-backup")
ok <- file.copy(fh, backup, overwrite = TRUE)
if (!isTRUE(ok)) {
  message(
    "patch_rcpp_function_h: could not create backup (likely read-only library): ",
    backup,
    " - skip patch"
  )
  quit(status = 0L)
}
wr_ok <- TRUE
tryCatch(
  writeLines(full_function_h, fh, useBytes = TRUE),
  error = function(e) {
    wr_ok <<- FALSE
    message(
      "patch_rcpp_function_h: could not write Function.h (likely read-only library): ",
      conditionMessage(e),
      " - restored backup and skipped patch"
    )
  }
)
if (!wr_ok) {
  try(file.copy(backup, fh, overwrite = TRUE), silent = TRUE)
  quit(status = 0L)
}

message(
  "patch_rcpp_function_h: replaced Function.h with full upstream conditional version. ",
  "expected branch=", exp_branch,
  if (is.na(det_branch)) "; detected branch=unknown" else paste0("; detected branch=", det_branch),
  "; backup=", backup
)
