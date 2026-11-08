//! PR/PR-BM homogeneous roots, native GERG-2008 T-density properties and (A7)
//! the NRTL excess model with a gamma–phi liquid/vapor package.
//! These are homogeneous evaluators; stability and flash live in `thermo-solvers`.
mod activity;
mod cubic;
mod nrtl;
pub use activity::*;
pub use cubic::*;
mod gerg;
pub use gerg::*;
pub use nrtl::*;
