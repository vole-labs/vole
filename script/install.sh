#!/usr/bin/env bash
#
# One-shot setup for vole on a fresh machine:
#   1. install system build dependencies (cmake, git, compiler, OpenSSL)
#   2. fetch the EMP toolkit git submodules (emp-tool / emp-ot 0.3.0, pinned
#      in .gitmodules; vole uses only emp-tool + emp-ot, no emp-zk)
#   3. build vole (this repo) against the submodules
#
# Supports Ubuntu/Debian (apt) and macOS (Homebrew).
#
# Usage:
#   script/install.sh
#
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

log()  { printf '\n\033[1;32m[install]\033[0m %s\n' "$*"; }
die()  { printf '\n\033[1;31m[install:error]\033[0m %s\n' "$*" >&2; exit 1; }
ncpu() { nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4; }

# ---- 0. platform + privilege ------------------------------------------------
OS="$(uname -s)"
SUDO=""
[ "$(id -u)" -eq 0 ] || SUDO="sudo"

OPENSSL_FLAG=()

# ---- 1. system dependencies -------------------------------------------------
case "$OS" in
  Linux)
    command -v apt-get >/dev/null 2>&1 || \
      die "Only apt-based Linux is scripted. Install cmake/git/g++/libssl-dev manually, then re-run with the EMP steps."
    log "Installing system packages via apt"
    $SUDO apt-get update -y
    $SUDO apt-get install -y build-essential cmake git libssl-dev
    ;;
  Darwin)
    command -v brew >/dev/null 2>&1 || \
      die "Homebrew not found. Install it from https://brew.sh and re-run."
    log "Installing system packages via Homebrew"
    brew install cmake git openssl@3 || true
    OPENSSL_FLAG=( -DOPENSSL_ROOT_DIR="$(brew --prefix openssl@3)" )
    ;;
  *)
    die "Unsupported OS '$OS' (supports Linux/apt and macOS/brew)."
    ;;
esac

# ---- 2. EMP toolkit (git submodules, pinned) --------------------------------
log "Fetching EMP submodules (emp-tool / emp-ot 0.3.0)"
git -C "$REPO_ROOT" submodule update --init

# ---- 3. build vole -------------------------------------------------------
log "Building vole"
cd "$REPO_ROOT"
rm -rf build
cmake -S . -B build "${OPENSSL_FLAG[@]}"
cmake --build build -j"$(ncpu)"

log "Done."
echo "  Correctness tests : (cd $REPO_ROOT/build && ctest)"
echo "  Benchmarks        : $REPO_ROOT/build/bin/bench_cvole 1 12345 z2k 4 22  (and party 2)"
echo "  Rebuild only      : cmake --build build -j"
