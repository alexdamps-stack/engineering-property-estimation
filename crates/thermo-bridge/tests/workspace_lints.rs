//! Guards the hand-copied `[lints]` table of this crate against drift from
//! `[workspace.lints]`.
//!
//! Cargo cannot combine `[lints] workspace = true` with a per-crate override,
//! and the workspace `unsafe_code = "forbid"` cannot be relaxed by an inner
//! `#[allow]`. The bridge therefore restates the workspace lints and changes
//! only `unsafe_code` to `deny` (the `ffi` module allows it). This test fails
//! when a workspace lint is added, removed or changed without the bridge
//! following.

use std::collections::BTreeMap;

/// `(table, key) -> value` for every `[<prefix><table>]` section of a TOML file.
fn lint_tables(toml: &str, prefix: &str) -> BTreeMap<(String, String), String> {
    let mut out = BTreeMap::new();
    let mut table: Option<String> = None;
    for line in toml.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('[').and_then(|l| l.strip_suffix(']')) {
            table = header.trim().strip_prefix(prefix).map(str::to_owned);
            continue;
        }
        if let (Some(t), Some((key, value))) = (&table, line.split_once('=')) {
            out.insert(
                (t.clone(), key.trim().to_owned()),
                value.trim().replace(' ', ""),
            );
        }
    }
    out
}

#[test]
fn bridge_lints_match_workspace_except_unsafe_code() {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../Cargo.toml");
    let own = concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml");
    let workspace = lint_tables(&std::fs::read_to_string(root).unwrap(), "workspace.lints.");
    let bridge = lint_tables(&std::fs::read_to_string(own).unwrap(), "lints.");
    assert!(!workspace.is_empty(), "no [workspace.lints.*] tables found");

    let unsafe_key = ("rust".to_owned(), "unsafe_code".to_owned());
    assert_eq!(
        workspace.get(&unsafe_key).map(String::as_str),
        Some("\"forbid\"")
    );
    assert_eq!(
        bridge.get(&unsafe_key).map(String::as_str),
        Some("\"deny\"")
    );

    let strip = |m: &BTreeMap<(String, String), String>| {
        let mut m = m.clone();
        m.remove(&unsafe_key);
        m
    };
    assert_eq!(
        strip(&bridge),
        strip(&workspace),
        "crates/thermo-bridge/Cargo.toml [lints] must restate [workspace.lints] (only unsafe_code differs)"
    );
}
