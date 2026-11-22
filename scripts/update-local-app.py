#!/usr/bin/env python3
"""Build, verify and atomically publish the local Linux/WSL desktop prototype."""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import uuid


ROOT = Path(__file__).resolve().parents[1]


def run(command, **kwargs):
    return subprocess.run(command, cwd=ROOT, check=True, **kwargs)


def git_command(*arguments):
    # Orca's Windows worktrees store a Windows Git directory in .git. Give
    # Linux Git an explicit translated path without rewriting shared Git files.
    command = ["git"]
    marker = ROOT / ".git"
    if marker.is_file():
        location = marker.read_text().strip().removeprefix("gitdir: ")
        if len(location) > 2 and location[1] == ":":
            location = run(["wslpath", "-u", location], capture_output=True,
                           text=True).stdout.strip()
            command += ["--git-dir=" + location, "--work-tree=" + str(ROOT)]
    return command + list(arguments)


def source_fingerprint():
    """Include dirty/untracked inputs; a Git commit alone misses local edits."""
    files = run(
        git_command("ls-files", "-z", "--cached", "--others", "--exclude-standard"),
        capture_output=True,
    ).stdout.split(b"\0")
    digest = hashlib.sha256()
    for name in sorted(set(files) - {b""}):
        path = ROOT / os.fsdecode(name)
        digest.update(len(name).to_bytes(8, "little"))
        digest.update(name)
        if path.is_file():
            digest.update(b"file\0")
            digest.update(hashlib.sha256(path.read_bytes()).digest())
        elif not path.exists():
            digest.update(b"deleted\0")
        else:
            raise RuntimeError(f"Unsupported source input: {path}")
    return digest.hexdigest()


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def is_current(app, fingerprint):
    try:
        manifest = json.loads((app / "current/build.json").read_text())
        return (
            manifest["source_root"] == str(ROOT)
            and manifest["source_sha256"] == fingerprint
            and manifest["binary_sha256"] == file_hash(app / "current/thermo-desktop")
        )
    except (OSError, ValueError, KeyError):
        return False


def write_notices(stage, metadata):
    notices = stage / "notices"
    notices.mkdir()
    shutil.copy2(ROOT / "packaging/windows/NOTICE.txt", notices)
    for name in ("CoolProp-LICENSE.txt", "Cantera-License.txt", "thermo-LICENSE.txt"):
        shutil.copy2(ROOT / "data/source-extracts" / name, notices)
    packages = {p["id"]: p for p in metadata["packages"]}
    nodes = {n["id"]: n for n in metadata["resolve"]["nodes"]}
    pending = [p["id"] for p in packages.values() if p["name"] == "thermo-desktop"]
    seen = set()
    while pending:
        item = pending.pop()
        if item in seen:
            continue
        seen.add(item)
        pending.extend(
            d["pkg"] for d in nodes[item]["deps"]
            if any(k["kind"] is None for k in d["dep_kinds"])
        )
    lines = ["Desktop prototype dependencies from Cargo.lock.",
             "Name Version | License | Repository", ""]
    for p in sorted((packages[i] for i in seen), key=lambda p: (p["name"], p["version"])):
        lines.append(f"{p['name']} {p['version']} | {p['license']} | {p['repository']}")
    (notices / "THIRD-PARTY-CRATES.txt").write_text("\n".join(lines) + "\n")


def update(app, fingerprint):
    cargo = shutil.which("cargo")
    if not cargo:
        cargo = str(Path.home() / ".cargo/bin/cargo")
        os.environ["PATH"] = str(Path(cargo).parent) + os.pathsep + os.environ.get("PATH", "")
    os.environ.setdefault("CARGO_TARGET_DIR", str(Path.home() / ".cache/thermo-target"))
    # Both commands run in ROOT, so rustup honors rust-toolchain.toml.
    run([cargo, "build", "-p", "thermo-desktop", "--bin", "thermo-desktop", "--locked"])
    rustc = str(Path(cargo).with_name("rustc")) if Path(cargo).is_absolute() else "rustc"
    host = next(line.split(": ", 1)[1] for line in
                run([rustc, "-vV"], capture_output=True, text=True).stdout.splitlines()
                if line.startswith("host: "))
    metadata = json.loads(run(
        [cargo, "metadata", "--format-version", "1", "--locked", "--filter-platform", host],
        capture_output=True, text=True,
    ).stdout)
    releases = app / "releases"
    releases.mkdir(exist_ok=True)
    current = app / "current"
    if current.exists() and not current.is_symlink():
        raise RuntimeError(f"Refusing to replace a non-symlink installation: {current}")
    with tempfile.TemporaryDirectory(prefix=".candidate-", dir=releases) as pending:
        stage = Path(pending)
        binary = stage / "thermo-desktop"
        shutil.copy2(Path(metadata["target_directory"]) / "debug/thermo-desktop", binary)
        if strip := shutil.which("strip"):
            run([strip, str(binary)])
        run([str(binary), "--self-test", str(stage / "self-test.json")])
        if json.loads((stage / "self-test.json").read_text()).get("status") != "passed":
            raise RuntimeError("Desktop acceptance self-test did not report a pass")
        write_notices(stage, metadata)
        shutil.copy2(ROOT / "LICENSE", stage / "LICENSE")
        shutil.copy2(ROOT / "packaging/windows/README.txt", stage / "README.txt")
        examples = stage / "examples"
        examples.mkdir()
        for example in (ROOT / "examples").glob("*.flowsheet.json"):
            shutil.copy2(example, examples)
        now = datetime.datetime.now(datetime.timezone.utc).isoformat()
        manifest = {
            "built_at": now,
            "source_root": str(ROOT),
            "source_commit": run(git_command("rev-parse", "HEAD"), capture_output=True,
                                 text=True).stdout.strip(),
            "source_sha256": fingerprint,
            "binary_sha256": file_hash(binary),
        }
        (stage / "build.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if source_fingerprint() != fingerprint:
            raise RuntimeError("Sources changed during the build; run the updater again")
        release = releases / uuid.uuid4().hex
        stage.rename(release)
        link = app / (".current-" + uuid.uuid4().hex)
        try:
            link.symlink_to(release.relative_to(app), target_is_directory=True)
            os.replace(link, current)
        finally:
            link.unlink(missing_ok=True)
    print(f"Updated desktop app: {current / 'thermo-desktop'}", flush=True)
    print("Open windows keep their existing build. Reopen the desktop shortcut to use this update.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Exit 1 if the installed app is stale")
    parser.add_argument("--app-dir", type=Path,
                        default=Path.home() / ".local/share/property-estimation-prototype")
    args = parser.parse_args()
    app = args.app_dir.expanduser().absolute()
    if args.check:
        current = is_current(app, source_fingerprint())
        print("Desktop app is current." if current else "Desktop app is missing or stale; run the updater.")
        return 0 if current else 1
    app.mkdir(parents=True, exist_ok=True)
    with (app / ".update.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        fingerprint = source_fingerprint()
        if is_current(app, fingerprint):
            print("Desktop app is already current.")
        else:
            update(app, fingerprint)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"Desktop update failed: {error}", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError) and error.stderr:
            print(error.stderr.decode() if isinstance(error.stderr, bytes) else error.stderr,
                  file=sys.stderr)
        sys.exit(2)
