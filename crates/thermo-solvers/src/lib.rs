//! Rachford–Rice, tangent-plane stability, pure-component saturation pressure and
//! two-phase TP flash for the PR/PR-BM cubic models.
//!
//! Stability is assessed from a finite, documented set of tangent-plane-distance
//! trials. Passing those trials is evidence, not a proof of global stability.
//! At most two phases (vapor–liquid or liquid–liquid by volume label) are solved;
//! an additional phase indicated after the split is reported as unsupported.
mod flash;
mod rachford_rice;
mod saturation;
mod stability;

pub use flash::*;
pub use rachford_rice::*;
pub use saturation::*;
pub use stability::*;
