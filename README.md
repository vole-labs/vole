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
| `mvole.h` | `MVoleFp<IO, FP, FPS>` | **n-party VOLE on primal LPN** (MVZK Protocol Π_nVOLE, after Le Mans): every pair runs a plain `VoleTriple`, a verifier reuses one programming seed towards all peers so its values coincide, the king regenerates them locally (`PrimalValueLocal`), and a final consistency check (fold with a fresh coin, zero-sharings, commit-and-open of `(u^i, Z^i_j)`, three equations) pins every verifier to a single value vector. Same correlation as `mcvole.h` at the primal rate, without a standalone commitment |

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
`test_ncvole_fp`, `test_mcvole_fp`, `test_mvole_fp`, `test_com_binding`, `test_mpfss_chi`,
`test_vole_f2k`, `test_vole_f2`, `bench_vole`, `bench_cvole`, `bench_mvole`,
`bench_mcvole`.

To build against installed emp-tool / emp-ot 1.0 packages instead of the
submodules, configure with `-DVOLE_USE_SYSTEM_EMP=ON`.

What this repo takes from EMP: `SoftSpoken<4>` (the only OT extension: the
correlated OTs for the GGM levels, `vole/base_cot.h`; `-DVOLE_COT_IKNP`
selects IKNP, `-DVOLE_COT_SOFTSPOKEN_K=k` another k), `CSW` (base OTs for
COPE), and emp-tool's
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

**n-party VOLE on primal LPN** (`mvole.h`; same process layout, king = `n`):

```sh
./build/bin/test_mvole_fp 0 12345 fp61 &
./build/bin/test_mvole_fp 1 12345 fp61 &
./build/bin/test_mvole_fp 2 12345 fp61 &
./build/bin/test_mvole_fp 3 12345 fp61        # king
```

Each verifier runs two extend rounds, sends a hash of its values to the king,
who compares them with its local regeneration, then all verifiers run the
Protocol-1 consistency check. A fifth argument `cheat` makes verifier 1 seed
its VOLE towards verifier 2 differently; every verifier must then abort (the
`mvole_fp_cheat` ctest case expects exit code 1). `bench_mvole <party> <port>
<field> <threads> <n_party> <rounds>` reports the per-round mesh, fold and
check times; `bench_mcvole` takes the same arguments (last one = log2 target)
for the committed dual-LPN `MCVoleFp`, so the two n-party paths compare
directly.

## References

- Chenkai Weng, Kang Yang, Jonathan Katz, and Xiao Wang. **Wolverine: Fast,
  Scalable, and Communication-Efficient Zero-Knowledge Proofs for Boolean and
  Arithmetic Circuits.** IEEE Symposium on Security and Privacy (S&P) 2021.
  IACR ePrint 2020/925. — the primal-LPN silent VOLE (`vole_triple.h`) and
  its parameters (`fp_default`), the base sVOLE and MPFSS consistency checks.
- Yunqing Sun, Hanlin Liu, Kang Yang, Yu Yu, Xiao Wang, and Chenkai Weng.
  **Committed Vector Oblivious Linear Evaluation and Its Applications.** 2025.
  — the committed VOLE (`cvole.h`, `com_matrix.h`): the dual-LPN expansion,
  the LPN-based commitment and its binding bound (Theorem 1), the multi-client
  and n-party variants (`ncvole.h`, `mcvole.h`).
