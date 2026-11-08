//! Untranslated PR and PR-BM homogeneous roots and molar caloric properties, and
//! (A7) the NRTL excess model with a gamma–phi liquid/vapor package.
//! These are homogeneous evaluators; stability and flash live in `thermo-solvers`.
mod activity;
mod cubic;
mod nrtl;
pub use activity::*;
pub use cubic::*;
pub use nrtl::*;
