//! In-memory project, its JSON file format, and solving.

use std::collections::BTreeMap;
use std::ops::ControlFlow;
use std::path::Path;

use flowsheet_core::{
    FeedCompositions, Flowsheet, PortRef, SolveProgress, Stream, StreamId, Unit, UnitId, UnitKind,
};
use property_packages::CubicPackage;
use serde::{Deserialize, Serialize};
use thermo_types::{
    COMPOSITION_SUM_TOLERANCE, ComponentId, Composition, CompositionBasis, NormalizationPolicy,
};

use crate::records::{KindRecord, SpecRecord};
use crate::{Fingerprints, PackageSelection, ProjectError, ResultSet};

/// Project file format version written by this program. Files declaring
/// any other version are rejected; there is no silent migration.
pub const FORMAT_VERSION: u64 = 1;

/// Feed mole fractions in package component order, per feed unit.
pub type FeedFractions = BTreeMap<UnitId, Vec<f64>>;

/// Display unit set stored with the project. Calculations are always SI.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DisplayUnitSet {
    /// °C, bar(a), kmol/h, kW.
    #[default]
    ProcessMetric,
    /// K, Pa(a), mol/s, W.
    Si,
}

/// A problem with project-level inputs that `flowsheet-core` does not hold,
/// such as feed compositions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputIssue {
    /// The affected unit.
    pub unit: UnitId,
    /// User-facing message naming the unit.
    pub message: String,
}

/// A solve that produced no result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolveFailure {
    /// The unit whose calculation failed, if the failure is attributable.
    pub unit: Option<UnitId>,
    /// User-facing cause.
    pub message: String,
    /// The solve was cancelled by the caller rather than failing.
    pub cancelled: bool,
}

/// Complete editable project: inputs, layout, display settings and the most
/// recent successful results.
#[derive(Debug, Clone, PartialEq)]
pub struct Project {
    /// Topology and unit specifications.
    pub flowsheet: Flowsheet,
    /// Property package selection for every stream.
    pub package: PackageSelection,
    /// Feed mole fractions in package order.
    pub feeds: FeedFractions,
    /// Canvas position (world coordinates) of each unit.
    pub layout: BTreeMap<UnitId, [f32; 2]>,
    /// Display unit set.
    pub display_units: DisplayUnitSet,
    /// Most recent successful results, if any. Check currency before use.
    pub results: Option<ResultSet>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PortRecord {
    unit: u64,
    port: usize,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct UnitRecord {
    id: u64,
    kind: KindRecord,
    tag: String,
    spec: SpecRecord,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StreamRecord {
    id: u64,
    name: String,
    from: PortRecord,
    to: PortRecord,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct FeedRecord {
    unit: u64,
    mole_fractions: Vec<f64>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LayoutRecord {
    unit: u64,
    x: f32,
    y: f32,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DisplayRecord {
    unit_set: DisplayUnitSet,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ProjectFile {
    format_version: u64,
    generator: String,
    package: PackageSelection,
    units: Vec<UnitRecord>,
    streams: Vec<StreamRecord>,
    feeds: Vec<FeedRecord>,
    layout: Vec<LayoutRecord>,
    display: DisplayRecord,
    results: Option<ResultSet>,
}

impl Project {
    /// An empty project with the given package.
    pub fn new(package: PackageSelection) -> Self {
        Self {
            flowsheet: Flowsheet::new(),
            package,
            feeds: FeedFractions::new(),
            layout: BTreeMap::new(),
            display_units: DisplayUnitSet::default(),
            results: None,
        }
    }

    /// Serializes to the versioned JSON project format. Results for units
    /// and streams that no longer exist are omitted.
    pub fn to_json(&self) -> Result<String, ProjectError> {
        let fs = &self.flowsheet;
        let port = |p: PortRef| PortRecord {
            unit: p.unit.0,
            port: p.port,
        };
        let results = self.results.clone().map(|mut r| {
            r.retain_existing(fs);
            r
        });
        let file = ProjectFile {
            format_version: FORMAT_VERSION,
            generator: concat!("thermo-desktop ", env!("CARGO_PKG_VERSION")).into(),
            package: self.package.clone(),
            units: fs
                .units()
                .map(|u| UnitRecord {
                    id: u.id.0,
                    kind: u.kind.into(),
                    tag: u.name.clone(),
                    spec: u.spec.into(),
                })
                .collect(),
            streams: fs
                .streams()
                .map(|s| StreamRecord {
                    id: s.id.0,
                    name: s.name.clone(),
                    from: port(s.from),
                    to: port(s.to),
                })
                .collect(),
            feeds: self
                .feeds
                .iter()
                .filter(|(id, _)| fs.unit(**id).is_some())
                .map(|(id, x)| FeedRecord {
                    unit: id.0,
                    mole_fractions: x.clone(),
                })
                .collect(),
            layout: self
                .layout
                .iter()
                .filter(|(id, _)| fs.unit(**id).is_some())
                .map(|(id, [x, y])| LayoutRecord {
                    unit: id.0,
                    x: *x,
                    y: *y,
                })
                .collect(),
            display: DisplayRecord {
                unit_set: self.display_units,
            },
            results,
        };
        let mut text =
            serde_json::to_string_pretty(&file).map_err(|e| ProjectError::Format(e.to_string()))?;
        text.push('\n');
        Ok(text)
    }

    /// Parses the versioned JSON project format. The version is checked
    /// before the content, so files from newer programs fail clearly.
    pub fn from_json(text: &str) -> Result<Self, ProjectError> {
        let value: serde_json::Value =
            serde_json::from_str(text).map_err(|e| ProjectError::Format(e.to_string()))?;
        match value
            .get("format_version")
            .and_then(serde_json::Value::as_u64)
        {
            Some(FORMAT_VERSION) => {}
            other => return Err(ProjectError::UnsupportedVersion(other)),
        }
        let file: ProjectFile =
            serde_json::from_value(value).map_err(|e| ProjectError::Format(e.to_string()))?;

        let mut units = Vec::new();
        for u in file.units {
            let spec = flowsheet_core::UnitSpec::try_from(u.spec)
                .map_err(|e| ProjectError::Invalid(format!("unit {}: {e}", u.tag)))?;
            units.push(Unit {
                id: UnitId(u.id),
                kind: u.kind.into(),
                name: u.tag,
                spec,
            });
        }
        let streams = file.streams.into_iter().map(|s| Stream {
            id: StreamId(s.id),
            name: s.name,
            from: PortRef::new(UnitId(s.from.unit), s.from.port),
            to: PortRef::new(UnitId(s.to.unit), s.to.port),
        });
        let flowsheet =
            Flowsheet::restore(units, streams).map_err(|e| ProjectError::Invalid(e.to_string()))?;

        let n = file.package.components.len();
        let mut feeds = FeedFractions::new();
        for f in file.feeds {
            let id = UnitId(f.unit);
            match flowsheet.unit(id) {
                Some(u) if u.kind == UnitKind::Feed => {}
                _ => {
                    return Err(ProjectError::Invalid(format!(
                        "composition given for {id}, which is not a feed"
                    )));
                }
            }
            if f.mole_fractions.len() != n {
                return Err(ProjectError::Invalid(format!(
                    "feed {id} has {} fractions for {n} components",
                    f.mole_fractions.len()
                )));
            }
            if feeds.insert(id, f.mole_fractions).is_some() {
                return Err(ProjectError::Invalid(format!("feed {id} appears twice")));
            }
        }
        let mut layout = BTreeMap::new();
        for l in file.layout {
            let id = UnitId(l.unit);
            if flowsheet.unit(id).is_none() || !l.x.is_finite() || !l.y.is_finite() {
                return Err(ProjectError::Invalid(format!(
                    "layout entry for {id} is invalid"
                )));
            }
            layout.insert(id, [l.x, l.y]);
        }
        Ok(Self {
            flowsheet,
            package: file.package,
            feeds,
            layout,
            display_units: file.display.unit_set,
            results: file.results,
        })
    }

    /// Writes the project file.
    pub fn save(&self, path: &Path) -> Result<(), ProjectError> {
        let text = self.to_json()?;
        std::fs::write(path, text).map_err(|e| ProjectError::Io(e.to_string()))
    }

    /// Reads a project file.
    pub fn load(path: &Path) -> Result<Self, ProjectError> {
        let text = std::fs::read_to_string(path).map_err(|e| ProjectError::Io(e.to_string()))?;
        Self::from_json(&text)
    }

    /// Current input fingerprints.
    pub fn fingerprints(&self) -> Fingerprints {
        Fingerprints::compute(&self.flowsheet, &self.feeds, &self.package)
    }

    /// Feed composition problems; the engine diagnostics cover everything else.
    pub fn input_issues(&self) -> Vec<InputIssue> {
        input_issues(&self.flowsheet, &self.feeds, &self.package)
    }

    /// Feed compositions with identities for the engine.
    pub fn feed_compositions(&self) -> Result<FeedCompositions, Vec<InputIssue>> {
        feed_compositions(&self.flowsheet, &self.feeds, &self.package)
    }
}

/// Feed composition problems of a flowsheet with the given package.
pub(crate) fn input_issues(
    fs: &Flowsheet,
    feeds: &FeedFractions,
    package: &PackageSelection,
) -> Vec<InputIssue> {
    let n = package.components.len();
    let mut out = Vec::new();
    for u in fs.units().filter(|u| u.kind == UnitKind::Feed) {
        let issue = |message: String| InputIssue {
            unit: u.id,
            message: format!("{}: {message}", u.name),
        };
        match feeds.get(&u.id) {
            None => out.push(issue("missing composition.".into())),
            Some(x) if x.len() != n => out.push(issue(format!(
                "composition has {} entries for {n} components.",
                x.len()
            ))),
            Some(x) if x.iter().any(|v| !v.is_finite() || *v < 0.0) => out.push(issue(
                "mole fractions must be finite and not negative.".into(),
            )),
            Some(x) => {
                let sum: f64 = x.iter().sum();
                if (sum - 1.0).abs() > COMPOSITION_SUM_TOLERANCE {
                    out.push(issue(format!(
                        "mole fractions sum to {sum:.10}; normalize them to 1."
                    )));
                }
            }
        }
    }
    out
}

/// Builds engine feed compositions, or returns the input issues.
pub(crate) fn feed_compositions(
    fs: &Flowsheet,
    feeds: &FeedFractions,
    package: &PackageSelection,
) -> Result<FeedCompositions, Vec<InputIssue>> {
    let issues = input_issues(fs, feeds, package);
    if !issues.is_empty() {
        return Err(issues);
    }
    let mut out = FeedCompositions::new();
    for u in fs.units().filter(|u| u.kind == UnitKind::Feed) {
        let fail = |e: &dyn std::fmt::Display| {
            vec![InputIssue {
                unit: u.id,
                message: format!("{}: {e}", u.name),
            }]
        };
        let entries = package
            .components
            .iter()
            .zip(&feeds[&u.id])
            .map(|(id, x)| Ok((ComponentId::new(id.clone()).map_err(|e| fail(&e))?, *x)))
            .collect::<Result<Vec<_>, Vec<InputIssue>>>()?;
        let c = Composition::new(
            CompositionBasis::Mole,
            entries,
            NormalizationPolicy::RequireUnitSum,
        )
        .map_err(|e| fail(&e))?;
        out.insert(u.id, c);
    }
    Ok(out)
}

/// Solves a project with an already built package (see
/// [`PackageSelection::build`]). `progress` is called before each unit;
/// returning [`ControlFlow::Break`] cancels between units.
///
/// The returned results carry the fingerprints of exactly the solved inputs.
/// A failure returns no results; the caller keeps any earlier results, whose
/// fingerprints show them stale where inputs changed.
pub fn solve_project(
    project: &Project,
    package: &CubicPackage,
    progress: impl FnMut(SolveProgress) -> ControlFlow<()>,
    solved_at_unix_s: u64,
) -> Result<ResultSet, SolveFailure> {
    let provenance = project.package.provenance().map_err(|e| SolveFailure {
        unit: None,
        message: e.to_string(),
        cancelled: false,
    })?;
    let feeds = project.feed_compositions().map_err(|issues| SolveFailure {
        unit: issues.first().map(|i| i.unit),
        message: issues
            .iter()
            .map(|i| i.message.as_str())
            .collect::<Vec<_>>()
            .join("\n"),
        cancelled: false,
    })?;
    let fingerprints = project.fingerprints();
    let calc = project
        .flowsheet
        .solve_with_progress(package, &feeds, progress)
        .map_err(|e| SolveFailure {
            unit: e.unit,
            message: e.message,
            cancelled: e.cancelled,
        })?;
    Ok(ResultSet::from_calculation(
        &calc,
        &project.flowsheet,
        &fingerprints,
        provenance,
        solved_at_unix_s,
    ))
}
