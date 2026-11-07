//! Dependency-free validated SI inputs and shared calculation vocabulary.
//!
//! Constructors reject invalid values. See `docs/calculation-contract.md` for
//! conventions and limits. No thermodynamic models or solvers are implemented.

mod composition;
mod error;
mod flow;
mod metadata;
mod quantity;

pub use composition::*;
pub use error::*;
pub use flow::*;
pub use metadata::*;
pub use quantity::*;
