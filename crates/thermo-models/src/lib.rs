//! PR/PR-BM homogeneous roots and native GERG-2008 T-density properties.
//! This is a homogeneous EOS evaluator; stability and flash live in `thermo-solvers`.
mod cubic;
pub use cubic::*;
mod gerg;
pub use gerg::*;
