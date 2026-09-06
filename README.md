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

Adding a field is a small change: write its header (a single-element type plus
its packed `value‖MAC` bundle), then include it in a test and add an
`else if (field == "…")` branch.

**Z_{2^k} (`z2k`).** The committed-VOLE stack is a linear map plus a single
multiply-by-Δ, so it is defined over any commutative ring — no field inverse is
needed. `z2k` instantiates the ring Z_{2^k} (SPDZ2k-style): the value
is in Z_{2⁶⁴}, the MAC is computed in Z_{2¹²⁸} (the extra s=64 bits are the
authentication-soundness margin), and Δ is confined to Z_{2⁶⁴}. All arithmetic
is native `unsigned __int128` (reduction mod 2¹²⁸ is hardware overflow); COPE's
bit-decomposition of Δ becomes the Gilboa shift-and-add product Δ·u = Σ aᵢ·2ⁱ.
A freshly sampled Δ is clamped to its 64-bit sub-ring by `restrict_delta()`
(a no-op for the field types). `z2k` is wired into the committed-VOLE tests
(`cvole_fp`, `ncvole_fp`, `mcvole_fp`); the plain primal-LPN `vole_fp` is
field-only.

> Soundness note: over a ring the universal-hash consistency checks (base sVOLE
> check, MPFSS batch check) rely on the s extra bits rather than field
> statistical soundness (the usual SPDZ2k heuristic). The commitment's binding
> argument (paper Theorem 1) is a field argument and does not cover the ring.

## Build

The EMP toolkit dependency (`emp-tool` and `emp-ot` at their **0.3.0**
releases) is vendored as git submodules under `thirdparty/`, pinned to the
exact commits the code was built against. Newer EMP releases (the 1.0 rewrite)
changed the OT API and do not build this code; `emp-zk` is not needed.

```sh
git clone --recursive <repo-url>        # or: git submodule update --init
cmake -S . -B build
cmake --build build -j
```

The only system requirements are a C++11 compiler, CMake ≥ 3.5, and OpenSSL
(`libssl-dev` / `brew install openssl@3`); `script/install.sh` installs those
on apt or Homebrew systems and then builds.

Binaries are placed in `build/bin/`: `test_vole_fp`, `test_cvole_fp`,
`test_ncvole_fp`, `test_mcvole_fp`, `test_com_binding`, `test_mpfss_chi`,
`bench_vole`, `bench_cvole`.

To build against an EMP 0.3.0 that is already installed instead of the
submodules, configure with `-DVOLE_USE_SYSTEM_EMP=ON` (add
`-DCMAKE_PREFIX_PATH=<prefix> -DCMAKE_FOLDER=<prefix>` if it is not in
`/usr/local`).

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
| fp61  | 1  | 2^20 | 2.10M | 5.69 s  | 2.74 s | 2.91 s | 2.71 |
| fp61  | 16 | 2^20 | 2.10M | 0.62 s  | 0.28 s | 0.30 s | 0.30 |
| f2k   | 1  | 2^20 | 2.10M | 7.46 s  | 4.11 s | 3.28 s | 3.56 |
| f2k   | 16 | 2^20 | 2.10M | 0.77 s  | 0.35 s | 0.34 s | 0.37 |
| fp107 | 1  | 2^20 | 2.10M | 12.69 s | 4.99 s | 7.61 s | 6.05 |
| fp107 | 16 | 2^20 | 2.10M | 1.15 s  | 0.42 s | 0.66 s | 0.55 |
| z2k   | 1  | 2^20 | 2.10M | 6.78 s  | 3.99 s | 2.72 s | 3.23 |
| z2k   | 16 | 2^20 | 2.10M | 0.75 s  | 0.39 s | 0.29 s | 0.36 |
| fp61  | 16 | 2^24 | 33.6M | 8.87 s  | 4.87 s | 3.49 s | 0.26 |
| f2k   | 16 | 2^24 | 33.6M | 11.44 s | 6.54 s | 3.96 s | 0.34 |

A target of 2^k needs two extend rounds (each round yields `n − M` usable
outputs), hence 2.10M correlations for 2^20. The commitment's cost is
dominated by per-column PRP generation and the `col_weight` multiply-adds, so
it is nearly independent of `n_com`.

The plain primal VOLE (`bench_vole fp61`, Wolverine parameters, 50M
correlations) runs at 0.028 / 0.011 / 0.0070 µs per correlation with 1 / 4 / 16
threads, level with emp-zk's `VoleTriple` single-threaded; with emp-ot's
`ferret_b13` parameter set (`bench_vole ... fp61 <threads> <batches> ferret_b13`)
it is 0.026 / 0.0090 / 0.0054 µs.
