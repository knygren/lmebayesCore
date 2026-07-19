#ifndef GLMBAYES_PACKAGE_NS_H
#define GLMBAYES_PACKAGE_NS_H

// R namespace for R callbacks that still live in glmbayesCore (EnvelopeOpt,
// EnvelopeSort, glmbfamfunc, rNormal_reg.wfit, rgamma_ct, rglmb),
// system.file("cl", ...), etc. These functions were removed from
// lmebayesCore's R/ (Stage 1 deduplication) and are now resolved from
// glmbayesCore's namespace instead.
// Override at compile time with -DGLMBAYES_R_NS=\"otherpack\" if needed.
#ifndef GLMBAYES_R_NS
#define GLMBAYES_R_NS "glmbayesCore"
#endif

#endif
