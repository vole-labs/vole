# Committed Vector Oblivious Linear Evaluation

A header-only, **field-templated** Vector Oblivious Linear Evaluation (VOLE)
library. The same code runs over many algebraic structures — prime fields, the
binary extension field GF(2¹²⁸), and the ring Z_{2^k} (SPDZ2k-style) — by
instantiating an element type, rather than being hardwired to a single field.

It implements the committed-VOLE (C-VOLE): a VOLE in
which the sender commits to its input vector and can then run VOLE with many
receivers, all provably using the same committed input.

## What's implemented (`vole/`)

| Header | Class | What it is |
|---|---|---|
| `vole_triple.h` | `VoleTriple<IO, FP, FPS>` | plain (primal-LPN) silent VOLE |
| `vole_f2k.h` | `RVole<IO, FP, FPS>`, `F2kVole<IO>`, `FerretVole<IO>` | Ferret-COT-style interface over `VoleTriple`: buffered random VOLE of any length, chosen-input VOLE with one-message derandomisation. `F2kVole` is the GF(2^128) instance with `block` arguments; `FerretVole` is the **F_2-value** instance (correlated OT, `R = K ⊕ b·Δ`, choice bit in the LSB) with Ferret's `rcot` / `send_cot` / `recv_cot` |
| `cvole.h` | `CVoleFp<IO, FP, FPS>` | **committed VOLE** (dual-LPN): VOLE on `x = H·A·e_u` plus the LPN commitment `com = H_u·e_u + H_r·e_r` (paper Sec. 4.1; `com_matrix.h`), with the consistency check `Hash(M[com]) == Hash(K[com] + com·Δ)` |
| `ncvole.h` | `ProgNCVoleFp<IO, FP, FPS>` | multi-client C-VOLE: one committer reuses a single committed input across many verifiers |
| `mcvole.h` | `MCVoleFp<IO, FP, FPS>` | **n-party pairwise committed VOLE**: a king seeds every party, reproduces and publishes each commitment locally, then every pair runs a committed VOLE|

## Fields (`vole/fields/`)

The tests support three fields, selected at run time:

| Name | Structure | Types `<FQ, FP, FPS>` |
|---|---|---|
| `fp61`  | F_p, p = 2⁶¹−1 (Mersenne) | `FP61, FP61, FP61x2` |
| `fp107` | F_p, p = 2¹⁰⁷−1 (Mersenne) | `FP107, FP107, FP107x2` |
| `f2k`   | GF(2¹²⁸) (binary extension) | `FP2x128, FP2x128, FP2x128x2` |
| `z2k`   | ring Z_{2^k}, k=64 (MAC in Z_{2¹²⁸}, Δ in Z_{2⁶⁴}) | `Z2k64, Z2k64, Z2k64x2` |
| `f2`    | values in F_2, keys/MACs/Δ in GF(2¹²⁸) (Ferret's correlated OT; one block per output, bit in the LSB) | `F2kKey, F2kKey, F2Auth` (`fields/f2.h`) |

Adding a field is a small change: write its header (a single-element type plus
its packed `value‖MAC` bundle), then include it in a test and add an
`else if (field == "…")` branch.

**Z_{2^k} (`z2k`).** The committed-VOLE stack is a linear map plus a single
multiply-by-Δ, so it is defined over any commutative ring — no field inverse is
needed. `z2k` instantiates the ring Z_{2^k} (SPDZ2k-style, k = s = 64): every
slot is a 128-bit `unsigned __int128` and all arithmetic is mod 2¹²⁸ (hardware
overflow); Δ is confined to Z_{2⁶⁴} (a freshly sampled Δ is clamped by
`restrict_delta()`, a no-op for the field types); COPE's bit-decomposition of
Δ becomes the Gilboa shift-and-add product Δ·u = Σ aᵢ·2ⁱ. `z2k` runs in every
driver and test (`vole_fp`, `cvole_fp`, `ncvole_fp`, `mcvole_fp`,
`com_binding`, `mpfss_chi`).

> **Value domain.** The value slot holds a full 128-bit residue, but only its
> low 64 bits are authenticated: with Δ ∈ Z_{2⁶⁴} the MAC relation
> M = K + x·Δ (mod 2¹²⁸) tolerates changes to the high 64 bits of x with
> noticeable probability. Consumers must treat `x mod 2⁶⁴` as the value and
> the high bits as an unauthenticated lift (as in SPDZ2k). Soundness of the
> checks is 2⁻⁶⁴ for the low 64 bits: the base-sVOLE check uses independent
> uniform ring coefficients (a polynomial hash degenerates in Z_{2^k}), the
> MPFSS checks use the MozZ2karella GGM-tag check plus binary coefficients.
>
> **Commitment.** The binding argument (paper Theorem 1) is a counting bound;
> over Z_{2^k} it reduces mod 2, so `n_com` is derived with one bit per row
> (about 10× the field value: 9536 rows at 2^20) and |e_r| is widened to
> 128·2⁸ so the dual-LPN hiding ratio n_com/N_com stays below the field's.

## Build

The EMP toolkit dependency is vendored as git submodules under `thirdparty/`:
`emp-tool` (upstream main, 1.0 line) and `emp-ot` (the `carlweng/emp-ot`
fork of the 1.0 line, which adds multithreaded silent sVOLE). Both build in
tree; no system-wide EMP install is needed. The code is C++20.

```sh
git clone --recursive <repo-url>        # or: git submodule update --init
cmake -S . -B build
cmake --build build -j
```

System requirements: a C++20 compiler (GCC 11+ / Clang 14+), CMake ≥ 3.21,
OpenSSL (`libssl-dev` / `brew install openssl@3`); `script/install.sh` installs
those on apt or Homebrew systems and then builds.

Binaries are placed in `build/bin/`: `test_vole_fp`, `test_cvole_fp`,
`test_ncvole_fp`, `test_mcvole_fp`, `test_com_binding`, `test_mpfss_chi`,
`test_vole_f2k`, `test_vole_f2`, `bench_vole`, `bench_cvole`.

To build against installed emp-tool / emp-ot 1.0 packages instead of the
submodules, configure with `-DVOLE_USE_SYSTEM_EMP=ON`.

What this repo takes from EMP: `SoftSpoken` (correlated OTs for the GGM
levels, `vole/base_cot.h`), `CSW` (base OTs for COPE), and emp-tool's
primitives (`PRG`, `PRP`, `CCRH`, `Hash`, `NetIO`, `ThreadPool`, `gfmul`). The
GGM node expander and the precomputed-OT layer that emp-ot 1.0 dropped are
vendored in `vole/twokeyprp.h` and `vole/preot.h`.

## Running the tests

Each test is a networked program; pick the field with the optional last
argument (`fp61` | `fp107` | `f2k` | `z2k`; defaults to `fp61`; `z2k` is
available for the committed-VOLE tests). Parties connect over localhost.

**Two-party (plain VOLE, committed VOLE):**

```sh
# party 1 (ALICE) and party 2 (BOB), same port and field
./build/bin/test_cvole_fp 1 12345 f2k &
./build/bin/test_cvole_fp 2 12345 f2k
```

`test_vole_fp` is invoked the same way.

**Multi-client committed VOLE** (server = party 0, then clients 1..n):

```sh
./build/bin/test_ncvole_fp 0 12345 fp107 &   # server (committer)
./build/bin/test_ncvole_fp 1 12345 fp107 &   # client 1 (verifier)
./build/bin/test_ncvole_fp 2 12345 fp107     # client 2 (verifier)
```

`extend_inplace_send/recv/local` return the number of usable VOLE outputs,
`x_usable(rounds) = rounds·(n − M)`; the last `M` entries of the x buffer are the
next round's base-sVOLE seed and are zeroed in the caller's buffer after being
copied internally, so they cannot be consumed as correlations.

Each test checks the VOLE relation `M = K + value·Δ`, the committed-VOLE
consistency check, and (for the committer) that the committed input is
reproducible / reused across sessions and clients. `test_cvole_fp` additionally
runs a local "king" with the committer's seed and asserts it reproduces the
committer's `x` / `com` values without interaction (the basis for `mcvole.h`).

**n-party committed VOLE** (parties `0..n`, king = `n`; the test fixes `n = 3`,
so four processes). All use the same port and field:

```sh
./build/bin/test_mcvole_fp 0 12345 f2k &   # VOLE party 0
./build/bin/test_mcvole_fp 1 12345 f2k &   # VOLE party 1
./build/bin/test_mcvole_fp 2 12345 f2k &   # VOLE party 2
./build/bin/test_mcvole_fp 3 12345 f2k     # P_king
```

The king publishes every party's commitment; each pairwise verifier checks its
key share against the published commitment. Success = every party completes
without a check failure (set `MCVOLE_DEBUG=1` to trace the per-peer schedule).

## Performance

`bench_cvole <party> <port> <field> <threads> <log2_target>` runs both parties
of the committed VOLE and splits the extend time into the VOLE pipeline
(COT + MPFSS(`e_u`) + `H`) and the commitment pipeline (COT + MPFSS(`e_r`) +
`[H_u | H_r]`). Numbers from one AWS instance (32 vCPU, both parties on
localhost, no network shaping), t = 224, rate 1/7, `col_weight` 32, `n_com`
derived per field.

| Field | Threads | Target | Correlations | Total | VOLE | Commitment | µs/corr |
|---|---|---|---|---|---|---|---|
| fp61  | 1  | 2^20 | 2.10M | 3.63 s  | 2.07 s | 1.53 s | 1.73 |
| fp61  | 16 | 2^20 | 2.10M | 0.42 s  | 0.24 s | 0.15 s | 0.20 |
| f2k   | 1  | 2^20 | 2.10M | 5.74 s  | 3.71 s | 1.98 s | 2.74 |
| f2k   | 16 | 2^20 | 2.10M | 0.56 s  | 0.29 s | 0.21 s | 0.27 |
| fp107 | 1  | 2^20 | 2.10M | 11.09 s | 4.47 s | 6.55 s | 5.29 |
| fp107 | 16 | 2^20 | 2.10M | 1.03 s  | 0.40 s | 0.57 s | 0.49 |
| z2k   | 1  | 2^20 | 2.10M | 6.13 s  | 3.55 s | 2.52 s | 2.93 |
| z2k   | 16 | 2^20 | 2.10M | 0.70 s  | 0.37 s | 0.26 s | 0.33 |
| fp61  | 16 | 2^24 | 33.6M | 6.92 s  | 4.59 s | 1.88 s | 0.21 |
| f2k   | 16 | 2^24 | 33.6M | 9.56 s  | 6.26 s | 2.45 s | 0.28 |

A target of 2^k needs two extend rounds (each round yields `n − M` usable
outputs), hence 2.10M correlations for 2^20. The commitment's cost is
dominated by per-column PRP generation and the `col_weight` multiply-adds, so
it is nearly independent of `n_com` for the fields; over `z2k` the dense `H_r`
block (32768 columns × 9536 rows) is what makes the commitment dearer.

The plain primal VOLE (`bench_vole`, 50M correlations, µs per correlation at
1 / 4 / 16 threads; buffers allocated outside the timed region):

| Instance | 1 thread | 4 threads | 16 threads | reference (measured against emp 0.3.0 on the same box) |
|---|---|---|---|---|
| `fp61` (Wolverine params) | 0.0174 | 0.0061 | 0.0038 | emp-zk `VoleTriple`: 0.025 / 0.0116 / 0.0059 (its LPN uses one extra worker) |
| `f2k` GF(2^128) values | 0.0245 | 0.0088 | 0.0052 | |
| `f2` F_2 values (Ferret COT, `fp_ferret_f2`) | 0.0180 | 0.0052 | 0.0026 | emp-ot 0.3.0 `FerretCOT`, same params: 0.0262 / 0.0073 / 0.0029 |
