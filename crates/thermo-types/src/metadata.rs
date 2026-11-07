/// Exact model selection; a family name alone is insufficient provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelMetadata {
    pub identifier: String,
    pub implementation_version: String,
    /// Alpha function, mixing rule, translation and other formulation choices.
    pub formulation: String,
}

/// Dataset revision and source locator; this record does not qualify a source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataMetadata {
    pub identifier: String,
    pub version: String,
    pub source: String,
}

/// Caloric zero/reference convention, distinct from gas-volume reference conditions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CaloricReference {
    pub identifier: String,
    /// Must specify h/s zeros, reference state, ideal contributions and units.
    pub definition: String,
}

/// No default: unassessed values must not imply a qualified domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidityStatus {
    WithinDeclaredDomain,
    Extrapolated,
    Unassessed,
}

/// Per-phase labels; no phase existence or stability is inferred by this enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    Vapor,
    Liquid,
    Supercritical,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiagnosticSeverity {
    Information,
    Warning,
}

/// Supplementary diagnostics never turn a failed calculation into success.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    pub code: String,
    pub severity: DiagnosticSeverity,
    pub message: String,
}
