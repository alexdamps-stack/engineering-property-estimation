//! Flowsheet topology: units, streams, and structural diagnostics.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::unit::{
    PortDirection, PortRef, PortSpec, SpecIssue, SpecProblem, Unit, UnitId, UnitKind, UnitSpec,
};

/// Stable identifier of a material stream within one [`Flowsheet`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StreamId(pub u64);

impl fmt::Display for StreamId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "S{}", self.0)
    }
}

/// A material stream from an outlet port to an inlet port.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Stream {
    /// Identifier, unique within its flowsheet.
    pub id: StreamId,
    /// User-facing stream name.
    pub name: String,
    /// Source outlet port.
    pub from: PortRef,
    /// Destination inlet port.
    pub to: PortRef,
}

/// Reasons a connection request is rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConnectError {
    /// A referenced unit does not exist.
    UnknownUnit(UnitId),
    /// A referenced port index does not exist on its unit.
    UnknownPort(PortRef),
    /// The source port is not an outlet.
    SourceNotOutlet(PortRef),
    /// The destination port is not an inlet.
    DestinationNotInlet(PortRef),
    /// The port already carries a stream.
    PortOccupied(PortRef),
    /// Source and destination belong to the same unit.
    SelfConnection(UnitId),
}

impl fmt::Display for ConnectError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConnectError::UnknownUnit(u) => write!(f, "unit {u} does not exist"),
            ConnectError::UnknownPort(p) => write!(f, "unit {} has no port {}", p.unit, p.port),
            ConnectError::SourceNotOutlet(_) => f.write_str("a stream must start at an outlet"),
            ConnectError::DestinationNotInlet(_) => f.write_str("a stream must end at an inlet"),
            ConnectError::PortOccupied(_) => f.write_str("the port is already connected"),
            ConnectError::SelfConnection(_) => f.write_str("a unit cannot be connected to itself"),
        }
    }
}

impl std::error::Error for ConnectError {}

/// Category of a structural diagnostic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiagnosticKind {
    /// A port has no stream attached.
    UnconnectedPort(PortRef),
    /// Required specification fields are unset.
    MissingSpecification(UnitId, Vec<&'static str>),
    /// A specified value is outside its admissible range.
    InvalidSpecification(UnitId, SpecIssue),
    /// The units form a recycle loop; recycles are planned for A10.
    RecycleUnsupported(Vec<UnitId>),
    /// The flowsheet contains no feed.
    NoFeed,
}

/// A structural problem that prevents calculation of the flowsheet.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    /// Machine-readable category.
    pub kind: DiagnosticKind,
    /// User-facing message.
    pub message: String,
}

/// Pre-calculation status of a unit (proposal section 4.6).
///
/// The calculated states (solved, warning, failed, stale) are added with the
/// flowsheet engine (A8); until then no unit is ever reported as solved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnitStatus {
    /// A specification is missing or invalid, or a port is unconnected.
    Underspecified,
    /// Complete and ready to be solved.
    Ready,
}

/// GUI-independent flowsheet topology and specifications.
#[derive(Debug, Clone, Default)]
pub struct Flowsheet {
    units: BTreeMap<UnitId, Unit>,
    streams: BTreeMap<StreamId, Stream>,
    next_unit: u64,
    next_stream: u64,
    tag_counters: BTreeMap<&'static str, u64>,
    revision: u64,
}

impl Flowsheet {
    /// Creates an empty flowsheet.
    pub fn new() -> Self {
        Self::default()
    }

    /// Monotonic counter incremented by every mutation.
    pub fn revision(&self) -> u64 {
        self.revision
    }

    /// Units in identifier order.
    pub fn units(&self) -> impl Iterator<Item = &Unit> {
        self.units.values()
    }

    /// Streams in identifier order.
    pub fn streams(&self) -> impl Iterator<Item = &Stream> {
        self.streams.values()
    }

    /// Looks up a unit.
    pub fn unit(&self, id: UnitId) -> Option<&Unit> {
        self.units.get(&id)
    }

    /// Looks up a stream.
    pub fn stream(&self, id: StreamId) -> Option<&Stream> {
        self.streams.get(&id)
    }

    /// Static description of a port, if it exists.
    pub fn port_spec(&self, port: PortRef) -> Option<PortSpec> {
        let unit = self.units.get(&port.unit)?;
        unit.kind.ports().get(port.port).copied()
    }

    /// The stream attached to a port, if any.
    pub fn stream_at(&self, port: PortRef) -> Option<&Stream> {
        self.streams
            .values()
            .find(|s| s.from == port || s.to == port)
    }

    /// Adds a unit with a generated tag and default specification.
    pub fn add_unit(&mut self, kind: UnitKind) -> UnitId {
        self.next_unit += 1;
        let id = UnitId(self.next_unit);
        let counter = self.tag_counters.entry(kind.tag_series()).or_insert(0);
        *counter += 1;
        let name = kind.tag(*counter);
        self.units.insert(
            id,
            Unit {
                id,
                kind,
                name,
                spec: kind.default_spec(),
            },
        );
        self.revision += 1;
        id
    }

    /// Removes a unit and every stream attached to it. Returns the removed unit.
    pub fn remove_unit(&mut self, id: UnitId) -> Option<Unit> {
        let unit = self.units.remove(&id)?;
        self.streams
            .retain(|_, s| s.from.unit != id && s.to.unit != id);
        self.revision += 1;
        Some(unit)
    }

    /// Renames a unit. Returns `false` if the unit does not exist.
    pub fn rename_unit(&mut self, id: UnitId, name: impl Into<String>) -> bool {
        let Some(unit) = self.units.get_mut(&id) else {
            return false;
        };
        unit.name = name.into();
        self.revision += 1;
        true
    }

    /// Replaces a unit's specification. Returns `false` if the unit does not
    /// exist or the specification variant does not belong to its kind.
    pub fn set_spec(&mut self, id: UnitId, spec: UnitSpec) -> bool {
        let Some(unit) = self.units.get_mut(&id) else {
            return false;
        };
        if std::mem::discriminant(&unit.kind.default_spec()) != std::mem::discriminant(&spec) {
            return false;
        }
        if unit.spec != spec {
            unit.spec = spec;
            self.revision += 1;
        }
        true
    }

    /// Connects an outlet port to an inlet port with a new stream.
    pub fn connect(&mut self, from: PortRef, to: PortRef) -> Result<StreamId, ConnectError> {
        let from_spec = self.checked_port(from)?;
        let to_spec = self.checked_port(to)?;
        if from_spec.direction != PortDirection::Outlet {
            return Err(ConnectError::SourceNotOutlet(from));
        }
        if to_spec.direction != PortDirection::Inlet {
            return Err(ConnectError::DestinationNotInlet(to));
        }
        if from.unit == to.unit {
            return Err(ConnectError::SelfConnection(from.unit));
        }
        for port in [from, to] {
            if self.stream_at(port).is_some() {
                return Err(ConnectError::PortOccupied(port));
            }
        }
        self.next_stream += 1;
        let id = StreamId(self.next_stream);
        self.streams.insert(
            id,
            Stream {
                id,
                name: id.to_string(),
                from,
                to,
            },
        );
        self.revision += 1;
        Ok(id)
    }

    /// Removes a stream. Returns the removed stream.
    pub fn disconnect(&mut self, id: StreamId) -> Option<Stream> {
        let stream = self.streams.remove(&id)?;
        self.revision += 1;
        Some(stream)
    }

    /// Renames a stream. Returns `false` if the stream does not exist.
    pub fn rename_stream(&mut self, id: StreamId, name: impl Into<String>) -> bool {
        let Some(stream) = self.streams.get_mut(&id) else {
            return false;
        };
        stream.name = name.into();
        self.revision += 1;
        true
    }

    /// Pre-calculation status of a unit, or `None` if it does not exist.
    pub fn unit_status(&self, id: UnitId) -> Option<UnitStatus> {
        let unit = self.units.get(&id)?;
        let open_port =
            (0..unit.kind.ports().len()).any(|i| self.stream_at(PortRef::new(id, i)).is_none());
        Some(if open_port || !unit.spec.issues().is_empty() {
            UnitStatus::Underspecified
        } else {
            UnitStatus::Ready
        })
    }

    fn checked_port(&self, port: PortRef) -> Result<PortSpec, ConnectError> {
        let unit = self
            .units
            .get(&port.unit)
            .ok_or(ConnectError::UnknownUnit(port.unit))?;
        unit.kind
            .ports()
            .get(port.port)
            .copied()
            .ok_or(ConnectError::UnknownPort(port))
    }

    /// Structural problems that must be resolved before calculation.
    ///
    /// An empty result means the topology and specifications are complete;
    /// it says nothing about whether a calculation would converge.
    pub fn diagnostics(&self) -> Vec<Diagnostic> {
        let mut out = Vec::new();
        if !self.units.values().any(|u| u.kind == UnitKind::Feed) {
            out.push(Diagnostic {
                kind: DiagnosticKind::NoFeed,
                message: "The flowsheet has no feed.".into(),
            });
        }
        for unit in self.units.values() {
            for (index, spec) in unit.kind.ports().iter().enumerate() {
                let port = PortRef::new(unit.id, index);
                if self.stream_at(port).is_none() {
                    out.push(Diagnostic {
                        kind: DiagnosticKind::UnconnectedPort(port),
                        message: format!("{}: port '{}' is not connected.", unit.name, spec.name),
                    });
                }
            }
            let (missing, invalid): (Vec<SpecIssue>, Vec<SpecIssue>) = unit
                .spec
                .issues()
                .into_iter()
                .partition(|i| i.problem == SpecProblem::Missing);
            if !missing.is_empty() {
                let fields: Vec<&'static str> = missing.iter().map(|i| i.field).collect();
                out.push(Diagnostic {
                    message: format!("{}: missing {}.", unit.name, fields.join(", ")),
                    kind: DiagnosticKind::MissingSpecification(unit.id, fields),
                });
            }
            for issue in invalid {
                out.push(Diagnostic {
                    message: format!("{}: {} {}.", unit.name, issue.field, issue.problem),
                    kind: DiagnosticKind::InvalidSpecification(unit.id, issue),
                });
            }
        }
        if let Some(cycle) = self.find_cycle() {
            let names: Vec<&str> = cycle
                .iter()
                .filter_map(|id| self.units.get(id).map(|u| u.name.as_str()))
                .collect();
            out.push(Diagnostic {
                message: format!(
                    "Recycle loop through {} — recycles are not supported yet (A10).",
                    names.join(" → ")
                ),
                kind: DiagnosticKind::RecycleUnsupported(cycle),
            });
        }
        out
    }

    /// Returns the units of one directed cycle, if any exists.
    fn find_cycle(&self) -> Option<Vec<UnitId>> {
        let mut successors: BTreeMap<UnitId, BTreeSet<UnitId>> = BTreeMap::new();
        for s in self.streams.values() {
            successors.entry(s.from.unit).or_default().insert(s.to.unit);
        }
        #[derive(Clone, Copy, PartialEq)]
        enum Mark {
            Active,
            Done,
        }
        let mut marks: BTreeMap<UnitId, Mark> = BTreeMap::new();
        for &start in self.units.keys() {
            if marks.contains_key(&start) {
                continue;
            }
            // Iterative depth-first search; `path` holds the active chain.
            let mut path: Vec<UnitId> = vec![start];
            let mut stack: Vec<Vec<UnitId>> = vec![Self::next_of(&successors, start)];
            marks.insert(start, Mark::Active);
            while let Some(pending) = stack.last_mut() {
                match pending.pop() {
                    Some(next) => match marks.get(&next) {
                        Some(Mark::Active) => {
                            let begin = path.iter().position(|&u| u == next).unwrap_or(0);
                            return Some(path[begin..].to_vec());
                        }
                        Some(Mark::Done) => {}
                        None => {
                            marks.insert(next, Mark::Active);
                            path.push(next);
                            stack.push(Self::next_of(&successors, next));
                        }
                    },
                    None => {
                        stack.pop();
                        if let Some(done) = path.pop() {
                            marks.insert(done, Mark::Done);
                        }
                    }
                }
            }
        }
        None
    }

    fn next_of(successors: &BTreeMap<UnitId, BTreeSet<UnitId>>, id: UnitId) -> Vec<UnitId> {
        successors
            .get(&id)
            .map(|s| s.iter().rev().copied().collect())
            .unwrap_or_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unit::{PressureSpec, ThermalSpec};

    fn out(unit: UnitId) -> PortRef {
        PortRef::new(unit, 0)
    }

    #[test]
    fn generated_tags_count_per_kind() {
        let mut fs = Flowsheet::new();
        let a = fs.add_unit(UnitKind::Feed);
        let b = fs.add_unit(UnitKind::Mixer);
        let c = fs.add_unit(UnitKind::Feed);
        let d = fs.add_unit(UnitKind::Heater);
        let e = fs.add_unit(UnitKind::Cooler);
        assert_eq!(fs.unit(a).unwrap().name, "Feed 1");
        assert_eq!(fs.unit(b).unwrap().name, "M-100");
        assert_eq!(fs.unit(c).unwrap().name, "Feed 2");
        // Heaters and coolers share the heat-exchanger series.
        assert_eq!(fs.unit(d).unwrap().name, "E-100");
        assert_eq!(fs.unit(e).unwrap().name, "E-101");
    }

    #[test]
    fn unit_status_requires_connections_and_valid_specs() {
        let mut fs = Flowsheet::new();
        let feed = fs.add_unit(UnitKind::Feed);
        let product = fs.add_unit(UnitKind::Product);
        assert_eq!(fs.unit_status(product), Some(UnitStatus::Underspecified));
        fs.connect(out(feed), PortRef::new(product, 0)).unwrap();
        assert_eq!(fs.unit_status(product), Some(UnitStatus::Ready));
        assert_eq!(fs.unit_status(feed), Some(UnitStatus::Underspecified));
        let valid = UnitSpec::Feed {
            temperature_k: Some(300.0),
            pressure_pa: Some(5.0e6),
            molar_flow_mol_s: Some(1.0),
        };
        fs.set_spec(feed, valid);
        assert_eq!(fs.unit_status(feed), Some(UnitStatus::Ready));
        fs.set_spec(
            feed,
            UnitSpec::Feed {
                temperature_k: Some(300.0),
                pressure_pa: Some(-1.0),
                molar_flow_mol_s: Some(1.0),
            },
        );
        assert_eq!(fs.unit_status(feed), Some(UnitStatus::Underspecified));
        assert!(fs.diagnostics().iter().any(|d| d.kind
            == DiagnosticKind::InvalidSpecification(
                feed,
                SpecIssue {
                    field: "pressure",
                    problem: SpecProblem::NotPositive
                }
            )));
        assert_eq!(fs.unit_status(UnitId(999)), None);
    }

    #[test]
    fn connect_enforces_direction_occupancy_and_self_loops() {
        let mut fs = Flowsheet::new();
        let feed = fs.add_unit(UnitKind::Feed);
        let heater = fs.add_unit(UnitKind::Heater);
        let heater_in = PortRef::new(heater, 0);
        let heater_out = PortRef::new(heater, 1);

        assert_eq!(
            fs.connect(heater_in, out(feed)),
            Err(ConnectError::SourceNotOutlet(heater_in))
        );
        assert_eq!(
            fs.connect(out(feed), heater_out),
            Err(ConnectError::DestinationNotInlet(heater_out))
        );
        assert_eq!(
            fs.connect(heater_out, heater_in),
            Err(ConnectError::SelfConnection(heater))
        );
        let bad = PortRef::new(feed, 7);
        assert_eq!(
            fs.connect(bad, heater_in),
            Err(ConnectError::UnknownPort(bad))
        );

        fs.connect(out(feed), heater_in).unwrap();
        let feed2 = fs.add_unit(UnitKind::Feed);
        assert_eq!(
            fs.connect(out(feed2), heater_in),
            Err(ConnectError::PortOccupied(heater_in))
        );
    }

    #[test]
    fn removing_a_unit_removes_attached_streams() {
        let mut fs = Flowsheet::new();
        let feed = fs.add_unit(UnitKind::Feed);
        let valve = fs.add_unit(UnitKind::Valve);
        let product = fs.add_unit(UnitKind::Product);
        fs.connect(out(feed), PortRef::new(valve, 0)).unwrap();
        fs.connect(PortRef::new(valve, 1), PortRef::new(product, 0))
            .unwrap();
        assert_eq!(fs.streams().count(), 2);
        fs.remove_unit(valve).unwrap();
        assert_eq!(fs.streams().count(), 0);
        assert_eq!(fs.units().count(), 2);
    }

    #[test]
    fn every_mutation_advances_the_revision() {
        let mut fs = Flowsheet::new();
        let r0 = fs.revision();
        let feed = fs.add_unit(UnitKind::Feed);
        let r1 = fs.revision();
        assert!(r1 > r0);
        let spec = UnitSpec::Feed {
            temperature_k: Some(300.0),
            pressure_pa: Some(5.0e6),
            molar_flow_mol_s: Some(1.0),
        };
        assert!(fs.set_spec(feed, spec));
        let r2 = fs.revision();
        assert!(r2 > r1);
        // Setting an identical spec is not a change.
        assert!(fs.set_spec(feed, spec));
        assert_eq!(fs.revision(), r2);
        assert!(fs.rename_unit(feed, "Natural gas"));
        assert!(fs.revision() > r2);
    }

    #[test]
    fn set_spec_rejects_variant_of_another_kind() {
        let mut fs = Flowsheet::new();
        let valve = fs.add_unit(UnitKind::Valve);
        let r = fs.revision();
        assert!(!fs.set_spec(
            valve,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::OutletTemperature(Some(280.0)),
                pressure_drop_pa: 0.0,
            }
        ));
        assert_eq!(fs.revision(), r);
    }

    #[test]
    fn complete_demonstration_topology_has_no_diagnostics() {
        // Mix two gas streams, cool, expand, separate (development plan §3).
        let mut fs = Flowsheet::new();
        let f1 = fs.add_unit(UnitKind::Feed);
        let f2 = fs.add_unit(UnitKind::Feed);
        let mix = fs.add_unit(UnitKind::Mixer);
        let cool = fs.add_unit(UnitKind::Cooler);
        let valve = fs.add_unit(UnitKind::Valve);
        let sep = fs.add_unit(UnitKind::Separator);
        let vap = fs.add_unit(UnitKind::Product);
        let liq = fs.add_unit(UnitKind::Product);
        for f in [f1, f2] {
            fs.set_spec(
                f,
                UnitSpec::Feed {
                    temperature_k: Some(300.0),
                    pressure_pa: Some(6.0e6),
                    molar_flow_mol_s: Some(1.0),
                },
            );
        }
        assert!(fs.set_spec(
            cool,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::OutletTemperature(Some(250.0)),
                pressure_drop_pa: 0.5e5,
            },
        ));
        assert!(fs.set_spec(
            valve,
            UnitSpec::Valve {
                pressure: PressureSpec::OutletPressure(Some(2.0e6)),
            },
        ));
        fs.connect(out(f1), PortRef::new(mix, 0)).unwrap();
        fs.connect(out(f2), PortRef::new(mix, 1)).unwrap();
        fs.connect(PortRef::new(mix, 2), PortRef::new(cool, 0))
            .unwrap();
        fs.connect(PortRef::new(cool, 1), PortRef::new(valve, 0))
            .unwrap();
        fs.connect(PortRef::new(valve, 1), PortRef::new(sep, 0))
            .unwrap();
        fs.connect(PortRef::new(sep, 1), PortRef::new(vap, 0))
            .unwrap();
        fs.connect(PortRef::new(sep, 2), PortRef::new(liq, 0))
            .unwrap();

        assert_eq!(fs.diagnostics(), Vec::new());
    }

    #[test]
    fn diagnostics_report_open_ports_missing_specs_and_no_feed() {
        let mut fs = Flowsheet::new();
        let valve = fs.add_unit(UnitKind::Valve);
        let kinds: Vec<_> = fs.diagnostics().into_iter().map(|d| d.kind).collect();
        assert!(kinds.contains(&DiagnosticKind::NoFeed));
        assert!(kinds.contains(&DiagnosticKind::UnconnectedPort(PortRef::new(valve, 0))));
        assert!(kinds.contains(&DiagnosticKind::UnconnectedPort(PortRef::new(valve, 1))));
        assert!(kinds.contains(&DiagnosticKind::MissingSpecification(
            valve,
            vec!["outlet pressure"]
        )));
    }

    #[test]
    fn recycle_loops_are_reported() {
        let mut fs = Flowsheet::new();
        let mix = fs.add_unit(UnitKind::Mixer);
        let heater = fs.add_unit(UnitKind::Heater);
        let sep = fs.add_unit(UnitKind::Separator);
        fs.connect(PortRef::new(mix, 2), PortRef::new(heater, 0))
            .unwrap();
        fs.connect(PortRef::new(heater, 1), PortRef::new(sep, 0))
            .unwrap();
        fs.connect(PortRef::new(sep, 2), PortRef::new(mix, 1))
            .unwrap();

        let cycle = fs
            .diagnostics()
            .into_iter()
            .find_map(|d| match d.kind {
                DiagnosticKind::RecycleUnsupported(c) => Some(c),
                _ => None,
            })
            .expect("recycle diagnostic");
        let mut sorted = cycle.clone();
        sorted.sort();
        assert_eq!(sorted, vec![mix, heater, sep]);
    }

    #[test]
    fn acyclic_branches_are_not_reported_as_recycles() {
        // Diamond: feed → separator → (vapor, liquid) → mixer.
        let mut fs = Flowsheet::new();
        let feed = fs.add_unit(UnitKind::Feed);
        let sep = fs.add_unit(UnitKind::Separator);
        let mix = fs.add_unit(UnitKind::Mixer);
        fs.connect(out(feed), PortRef::new(sep, 0)).unwrap();
        fs.connect(PortRef::new(sep, 1), PortRef::new(mix, 0))
            .unwrap();
        fs.connect(PortRef::new(sep, 2), PortRef::new(mix, 1))
            .unwrap();
        assert!(
            !fs.diagnostics()
                .iter()
                .any(|d| matches!(d.kind, DiagnosticKind::RecycleUnsupported(_)))
        );
    }
}
