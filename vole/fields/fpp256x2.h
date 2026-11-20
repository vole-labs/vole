#ifndef VOLE_FIELDS_FPP256X2_H__
#define VOLE_FIELDS_FPP256X2_H__
// FPP256x2: (value, MAC) pair over the NIST P-256 base field (64 bytes).
#include "vole/fields/fpp256.h"
#include "vole/fields/field_pair.h"
using FPP256x2 = FieldPair<FPP256>;
static_assert(sizeof(FPP256x2) == 64, "FPP256x2 must be 64 bytes");
#endif  // VOLE_FIELDS_FPP256X2_H__
