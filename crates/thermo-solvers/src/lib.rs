//! Rachford–Rice, tangent-plane stability, bounded pure/mixture bubble/dew boundaries and
//! two-phase TP/PH/PS flash for the PR/PR-BM cubic models.
//! Native GERG adds bracketed density inversion and dry-gas TPD/TP/PH integration
//! with an explicit exterior-root search and native caloric reference.
//! (A7) adds gamma–phi bubble/dew points and TP/PH flash for the NRTL activity package.
//!
//! Stability is assessed from a finite, documented set of tangent-plane-distance
//! trials. Passing those trials is evidence, not a proof of global stability.
//! The two-phase flashes solve at most two phases (vapor–liquid or liquid–liquid by
//! volume label) and report an additional phase indicated after the split as
//! unsupported; `multiphase` adds VLLE/LLE stability and TP/PH flash (#224).
mod activity;
mod equilibrium_model;
mod flash;
mod gerg_density;
mod gerg_flash;
mod gerg_ph;
mod mixture_saturation;
mod multiphase;
mod ph;
mod ps;
mod rachford_rice;
mod saturation;
mod stability;

pub use activity::*;
pub use flash::*;
pub use gerg_density::*;
pub use gerg_flash::*;
pub use gerg_ph::*;
pub use mixture_saturation::*;
pub use multiphase::*;
pub use ph::*;
pub use ps::*;
pub use rachford_rice::*;
pub use saturation::*;
pub use stability::*;
