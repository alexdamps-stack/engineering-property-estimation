//! Command-line behaviour of the `thermo-desktop` binary that needs no window.

use std::process::Command;

#[test]
fn version_to_an_unwritable_path_fails() {
    let dir = std::env::temp_dir().join(format!("thermo-desktop-cli-{}", std::process::id()));
    let missing = dir.join("no-such-directory").join("version.txt");
    let status = Command::new(env!("CARGO_BIN_EXE_thermo-desktop"))
        .arg("--version")
        .arg(&missing)
        .status()
        .expect("binary runs");
    assert!(!status.success(), "a lost version file must not exit 0");
}

#[test]
fn version_to_a_file_succeeds() {
    let dir = std::env::temp_dir().join(format!("thermo-desktop-cli-ok-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("version.txt");
    let status = Command::new(env!("CARGO_BIN_EXE_thermo-desktop"))
        .arg("--version")
        .arg(&path)
        .status()
        .expect("binary runs");
    assert!(status.success());
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.starts_with("thermo-desktop "), "{text}");
    let _ = std::fs::remove_dir_all(&dir);
}
