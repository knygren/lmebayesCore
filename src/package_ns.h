#ifndef GLMBAYES_PACKAGE_NS_H
#define GLMBAYES_PACKAGE_NS_H

// R namespace for package-local R callbacks (EnvelopeOpt, EnvelopeSort,
// glmbfamfunc, rNormal_reg.wfit, rgamma_ct), system.file("cl", ...), etc.
// Override at compile time with -DGLMBAYES_R_NS=\"otherpack\" if needed.
#ifndef GLMBAYES_R_NS
#define GLMBAYES_R_NS "glmbayesCore"
#endif

#endif
