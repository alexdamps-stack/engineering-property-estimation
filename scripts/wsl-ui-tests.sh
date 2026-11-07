#!/usr/bin/env bash
# Runs the headless desktop UI tests inside WSL (Linux).
#
# Use this on Windows machines where Smart App Control or another code
# integrity policy blocks freshly built, unsigned test executables. Linux
# binaries in WSL are not subject to that policy. Off-screen rendering for
# the snapshot tests uses Mesa's lavapipe Vulkan driver
# (`sudo apt install mesa-vulkan-drivers`).
#
# From Windows PowerShell, in the repository root (`--cd` needs an absolute path):
#   wsl -d Ubuntu --cd "$PWD" -- bash scripts/wsl-ui-tests.sh
# Extra arguments are passed to the test binary, for example a name filter:
#   wsl -d Ubuntu --cd "$PWD" -- bash scripts/wsl-ui-tests.sh snapshot
# Record new snapshot baselines after an intended visual change:
#   wsl -d Ubuntu --cd "$PWD" -- env UPDATE_SNAPSHOTS=1 bash scripts/wsl-ui-tests.sh
set -euo pipefail

# `wsl -- bash script` starts a non-login shell without rustup on PATH.
if ! command -v cargo >/dev/null 2>&1 && [ -f "$HOME/.cargo/env" ]; then
  # shellcheck disable=SC1091
  . "$HOME/.cargo/env"
fi

# Keep Linux build output out of the Windows `target/` directory and on the
# faster Linux file system.
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$HOME/.cache/thermo-target}"
# Force headless rendering through lavapipe so results do not depend on the
# host GPU passthrough.
export WGPU_BACKEND="${WGPU_BACKEND:-vulkan}"

cargo test -p thermo-desktop --test ui --locked -- "$@"
