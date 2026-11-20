#ifndef VOLE_FIELDS_FP25519X2_H__
#define VOLE_FIELDS_FP25519X2_H__
// FP25519x2: (value, MAC) pair over F_p, p = 2^255 - 19 (64 bytes).
#include "vole/fields/fp25519.h"
#include "vole/fields/field_pair.h"
using FP25519x2 = FieldPair<FP25519>;
static_assert(sizeof(FP25519x2) == 64, "FP25519x2 must be 64 bytes");
#endif  // VOLE_FIELDS_FP25519X2_H__
