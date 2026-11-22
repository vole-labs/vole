#ifndef VOLE_FIELDS_FP128X2_H__
#define VOLE_FIELDS_FP128X2_H__
// FP128x2: (value, MAC) pair over F_p, p = 2^128 - 159 (32 bytes).
#include "vole/fields/fp128.h"
#include "vole/fields/field_pair.h"
using FP128x2 = FieldPair<FP128>;
static_assert(sizeof(FP128x2) == 32, "FP128x2 must be 32 bytes");
#endif  // VOLE_FIELDS_FP128X2_H__
