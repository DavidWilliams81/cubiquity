#ifndef CUBIQUITY_VIEW_H
#define CUBIQUITY_VIEW_H

// Sanity check that build system is setting things up correctly.
#ifndef CUBIQUITY_APP_ENABLE_VIEW
#error "CUBIQUITY_APP_ENABLE_VIEW should be set when building view command"
#endif

#include "flags.h"

bool viewVolume(const flags::args& args);

#endif // CUBIQUITY_VIEW_H

