#!/bin/sh
# g++ wrapper that pins the code-generation target of the NIST GERG reference
# builds (prepare-gerg.py and the GERG density/equilibrium generators).
#
# The committed GERG fixtures were produced by a g++ whose default target is
# the generic x86-64 baseline. Some toolchains (for example the Ubuntu 26.04
# amd64v3 variant used by GitHub's ubuntu-latest runners, configured with
# --with-arch-64=x86-64-v3) default to AVX2/FMA; because GNU C++ contracts
# a*b+c into FMA by default, the reference values then differ in the last bits
# and the byte-exact reproduction fails. These flags restore the baseline
# target without editing the generators, whose hashes are recorded in the
# fixtures' provenance.
#
# Usage: put a directory containing a `g++` link to this script first on PATH,
# or call it directly. REAL_GXX selects the compiler (default /usr/bin/g++).
exec "${REAL_GXX:-/usr/bin/g++}" -march=x86-64 -mtune=generic -ffp-contract=off "$@"
