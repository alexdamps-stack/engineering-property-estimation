//! Stream transport properties (#225) in solved, saved and reloaded projects.
use std::ops::ControlFlow;

use flowsheet_project::{
    CubicVariantChoice, ENGINE_REVISION, Project, demonstration, solve_project,
};

#[test]
fn engine_revision_names_the_transport_versions() {
    assert!(ENGINE_REVISION.contains(thermo_models::TRANSPORT_MODEL_VERSION));
    assert!(ENGINE_REVISION.contains(thermo_data::TRANSPORT_DATA_VERSION));
}

#[test]
fn demonstration_streams_carry_transport_properties_through_save_and_reload() {
    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    let package = project.package.prepare().unwrap();
    let results = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
    let mut gas_phases = 0;
    let mut liquid_phases = 0;
    for stream in results.streams.values() {
        let Some(state) = &stream.state else { continue };
        for phase in &state.phases {
            let transport = phase.transport.as_ref().expect("evaluated");
            assert_eq!(transport.validity.as_deref(), Some("unassessed"));
            if phase.label == "liquid" {
                liquid_phases += 1;
                // The separator liquid holds methane far above 0.75 Tc.
                assert_eq!(transport.regime.as_deref(), Some("liquid"));
                assert_eq!(phase.dynamic_viscosity_pa_s, None);
                assert_eq!(phase.surface_tension_n_m, None);
                assert_eq!(transport.unavailable.len(), 3);
                assert!(transport.unavailable.iter().all(|r| r.contains("methane")));
            } else {
                gas_phases += 1;
                let mu = phase.dynamic_viscosity_pa_s.unwrap();
                let k = phase.thermal_conductivity_w_m_k.unwrap();
                // Natural-gas magnitudes at 240..350 K, 1.5..5 MPa.
                assert!((8e-6..3e-5).contains(&mu), "{mu}");
                assert!((0.015..0.06).contains(&k), "{k}");
                assert_eq!(phase.surface_tension_n_m, None);
                assert!(transport.unavailable.is_empty());
                assert_eq!(transport.method.as_deref(), Some("chung-1988 transport-v1"));
            }
        }
    }
    assert!(
        gas_phases >= 6 && liquid_phases >= 1,
        "{gas_phases} {liquid_phases}"
    );
    project.results = Some(results);
    let json = project.to_json().unwrap();
    assert!(json.contains("\"transport\""));
    let loaded = Project::from_json(&json).unwrap();
    assert_eq!(loaded.results, project.results);
    assert_eq!(loaded.to_json().unwrap(), json);
}
