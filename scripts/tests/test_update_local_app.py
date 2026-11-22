"""Exercise publication and failure recovery without compiling Rust or changing the real app."""
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


UPDATER = Path(__file__).resolve().parents[1] / "update-local-app.py"


class DesktopUpdateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="desktop-updater-test-")
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.repo, self.app, tools = root / "repo", root / "app", root / "tools"
        (self.repo / "scripts").mkdir(parents=True)
        tools.mkdir()
        shutil.copy2(UPDATER, self.repo / "scripts/update-local-app.py")
        for name in ("LICENSE", "Cargo.toml", "crates/example.rs", "packaging/windows/NOTICE.txt",
                     "packaging/windows/README.txt", "data/source-extracts/CoolProp-LICENSE.txt",
                     "data/source-extracts/Cantera-License.txt", "data/source-extracts/thermo-LICENSE.txt"):
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("test input\n")
        (self.repo / "examples").mkdir()
        self.env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"],
                        CARGO_TARGET_DIR=str(root / "build"))
        for key in ("UPDATER_TEST_BUILD_FAILURE", "UPDATER_TEST_SELF_TEST_FAILURE", "UPDATER_TEST_EDIT"):
            self.env.pop(key, None)
        self.tool(tools / "cargo", '''#!/usr/bin/env python3
import json, os, pathlib, sys
target = pathlib.Path(os.environ['CARGO_TARGET_DIR'])
if sys.argv[1] == 'build':
    if os.environ.get('UPDATER_TEST_BUILD_FAILURE'): sys.exit(23)
    if edit := os.environ.get('UPDATER_TEST_EDIT'): pathlib.Path(edit).write_text('changed during build')
    binary = target / 'debug/thermo-desktop'
    binary.parent.mkdir(parents=True, exist_ok=True)
    binary.write_text("#!/usr/bin/env python3\\nimport json,os,pathlib,sys\\nif os.environ.get('UPDATER_TEST_SELF_TEST_FAILURE'): sys.exit(24)\\npathlib.Path(sys.argv[2]).write_text(json.dumps({'status':'passed'}))\\n")
    binary.chmod(0o755)
elif sys.argv[1] == 'metadata':
    print(json.dumps({'target_directory':str(target),'packages':[{'id':'app','name':'thermo-desktop','version':'test','license':'MIT','repository':None}],'resolve':{'nodes':[{'id':'app','deps':[]}]}}))
''')
        self.tool(tools / "rustc", "#!/bin/sh\necho 'host: x86_64-unknown-linux-gnu'\n")
        self.tool(tools / "strip", "#!/bin/sh\nexit 0\n")
        for command in (["git", "init", "-q"], ["git", "add", "."],
                        ["git", "-c", "user.name=Updater Test", "-c", "user.email=test@example.invalid",
                         "commit", "-qm", "test inputs"]):
            subprocess.run(command, cwd=self.repo, env=self.env, check=True, capture_output=True)

    @staticmethod
    def tool(path, content):
        path.write_text(content)
        path.chmod(0o755)

    def invoke(self, check=False):
        command = ["python3", str(self.repo / "scripts/update-local-app.py"), "--app-dir", str(self.app)]
        if check:
            command.append("--check")
        return subprocess.run(command, cwd=self.repo, env=self.env, capture_output=True, text=True)

    def publish(self):
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.invoke(check=True).returncode, 0)
        return (self.app / "current").resolve(), (self.app / "current/build.json").read_bytes()

    def assert_preserved(self, previous):
        self.assertEqual((self.app / "current").resolve(), previous[0])
        self.assertEqual((self.app / "current/build.json").read_bytes(), previous[1])
        self.assertFalse(list((self.app / "releases").glob(".candidate-*")))

    def test_success_records_verified_binary_and_skips_unchanged_build(self):
        self.publish()
        manifest = json.loads((self.app / "current/build.json").read_text())
        self.assertEqual(manifest["source_root"], str(self.repo))
        self.assertEqual(json.loads((self.app / "current/self-test.json").read_text())["status"], "passed")
        self.env["UPDATER_TEST_BUILD_FAILURE"] = "1"
        self.assertEqual(self.invoke().returncode, 0)

    def test_untracked_edits_deletions_and_binary_tampering_are_detected(self):
        self.publish()
        extra = self.repo / "new-source.rs"
        extra.write_text("new input")
        self.assertEqual(self.invoke(check=True).returncode, 1)
        extra.unlink()
        self.assertEqual(self.invoke(check=True).returncode, 0)
        (self.repo / "crates/example.rs").unlink()
        self.assertEqual(self.invoke(check=True).returncode, 1)
        (self.repo / "crates/example.rs").write_text("test input\n")
        with (self.app / "current/thermo-desktop").open("ab") as binary:
            binary.write(b"changed")
        self.assertEqual(self.invoke(check=True).returncode, 1)

    def test_failed_build_keeps_previous_installation(self):
        previous = self.publish()
        (self.repo / "crates/example.rs").write_text("edited")
        self.env["UPDATER_TEST_BUILD_FAILURE"] = "1"
        self.assertEqual(self.invoke().returncode, 2)
        self.assert_preserved(previous)

    def test_failed_self_test_keeps_previous_installation(self):
        previous = self.publish()
        (self.repo / "crates/example.rs").write_text("edited")
        self.env["UPDATER_TEST_SELF_TEST_FAILURE"] = "1"
        self.assertEqual(self.invoke().returncode, 2)
        self.assert_preserved(previous)

    def test_sources_changing_during_build_are_not_published(self):
        previous = self.publish()
        source = self.repo / "crates/example.rs"
        source.write_text("edited")
        self.env["UPDATER_TEST_EDIT"] = str(source)
        result = self.invoke()
        self.assertEqual(result.returncode, 2)
        self.assertIn("Sources changed during the build", result.stderr)
        self.assert_preserved(previous)

    def test_competing_update_is_refused_without_replacing_the_app(self):
        previous = self.publish()
        (self.repo / "crates/example.rs").write_text("edited")
        with (self.app / ".update.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            self.assertEqual(self.invoke().returncode, 2)
        self.assert_preserved(previous)


if __name__ == "__main__":
    unittest.main()
