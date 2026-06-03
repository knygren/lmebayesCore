#ifndef GLMBAYES_PACKAGE_NS_H
#define GLMBAYES_PACKAGE_NS_H

// R namespace for glmbfamfunc, EnvelopeSort, system.file("cl", ...), etc.
// Override at compile time with -DGLMBAYES_R_NS=\"otherpack\" if needed.
#ifndef GLMBAYES_R_NS
#define GLMBAYES_R_NS "glmbayesCore"
#endif

#endif
