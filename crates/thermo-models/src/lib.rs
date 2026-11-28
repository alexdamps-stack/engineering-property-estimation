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

/// NASA7 caloric zero shared by the cubic and gamma–phi models.
pub const IDEAL_GAS_CALORIC_REFERENCE_ID: &str = "ideal-pure-zero-298.15K-1bar-cp-integrated-v1";

fn ideal_gas_caloric_reference() -> thermo_types::CaloricReference {
    thermo_types::CaloricReference {
        identifier: IDEAL_GAS_CALORIC_REFERENCE_ID.into(),
        definition: "NASA7 cp integrals; pure ideal h=s=0 at 298.15 K, 100000 Pa; ideal mixing entropy; J/mol and J/(mol K)".into(),
    }
}
mod eos_cg;
pub use eos_cg::{EosCg2021, EosCg2021State};
mod liquid_density;
pub use liquid_density::*;
mod transport;
pub use transport::*;
mod diffusivity;
pub use diffusivity::*;
