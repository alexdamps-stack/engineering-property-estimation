//! Reproducible homogeneous property calculations, with explicit failure output.

use std::{env, fs, process::ExitCode};

fn run() -> Result<(), String> {
    let args: Vec<_> = env::args().skip(1).collect();
    match args
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .as_slice()
    {
        [] | ["--help"] | ["-h"] => println!(
            "thermo-cli\n  components\n  calculate --input <request.json>\n\nHomogeneous PR/PR-BM roots only; no stability assessment or flash.\nSee docs/calculation-contract.md and examples/."
        ),
        ["components"] => {
            let catalog = thermo_data::Catalog::bundled().map_err(|e| e.to_string())?;
            println!(
                "{}",
                serde_json::to_string_pretty(&serde_json::json!({
                    "data_version": catalog.version,
                    "components": catalog.component_ids().collect::<Vec<_>>(),
                }))
                .map_err(|e| e.to_string())?
            );
        }
        ["calculate", "--input", path] => {
            let text = fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
            let request =
                serde_json::from_str(&text).map_err(|e| format!("invalid request JSON: {e}"))?;
            let result = property_packages::calculate(request).map_err(|e| e.to_string())?;
            println!(
                "{}",
                serde_json::to_string_pretty(&result).map_err(|e| e.to_string())?
            );
        }
        _ => return Err("unknown arguments; use --help".into()),
    }
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{}", serde_json::json!({"error": error}));
            ExitCode::FAILURE
        }
    }
}
