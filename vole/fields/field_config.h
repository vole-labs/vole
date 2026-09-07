#ifndef MPC_FIELD_CONFIG_H__
#define MPC_FIELD_CONFIG_H__

// Aggregates the supported field types plus the shared hashing / inner-product
// helpers (utils.h) for the templated VOLE components. Field selection is a
// runtime choice in the tests, which instantiate FP61 / FP107 / FP2x128
// directly; this header just makes those types and the helpers available.

#include "vole/fields/utils.h"
#include "vole/fields/fp61.h"
#include "vole/fields/fp61x2.h"
#include "vole/fields/fp107.h"
#include "vole/fields/fp107x2.h"
#include "vole/fields/fp2x128.h"
#include "vole/fields/z2k.h"
#include "vole/fields/f2.h"

#endif // MPC_FIELD_CONFIG_H__
