//! Headless A10 demonstration; optionally writes an input-only desktop project.
use flowsheet_project::{CubicVariantChoice, recycle_demonstration, solve_project};
use std::ops::ControlFlow;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let project = recycle_demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    if let Some(path) = std::env::args().nth(1) {
        project.save(std::path::Path::new(&path))?;
    }
    let package = project.package.prepare()?;
    let results = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0)
        .map_err(|e| e.message)?;
    println!("{}", serde_json::to_string_pretty(&results)?);
    Ok(())
}
