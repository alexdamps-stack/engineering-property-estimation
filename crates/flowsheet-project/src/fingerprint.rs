//! Input fingerprints for precise result invalidation.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use flowsheet_core::{Flowsheet, PortDirection, PortRef, StreamId, UnitId};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::records::{KindRecord, SpecRecord};
use crate::{FeedFractions, PackageSelection};

/// 64-bit FNV-1a digest of a canonical input description. It detects
/// changes; it is not a cryptographic hash. Stored as 16 hex digits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Fingerprint(pub u64);

impl fmt::Display for Fingerprint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:016x}", self.0)
    }
}

impl Serialize for Fingerprint {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Fingerprint {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let text = String::deserialize(d)?;
        if text.len() != 16 {
            return Err(serde::de::Error::custom(
                "fingerprint must have 16 hex digits",
            ));
        }
        u64::from_str_radix(&text, 16)
            .map(Fingerprint)
            .map_err(serde::de::Error::custom)
    }
}

#[derive(Clone, Copy)]
struct Fnv(u64);

impl Fnv {
    fn new() -> Self {
        Fnv(0xcbf2_9ce4_8422_2325)
    }
    fn bytes(&mut self, bytes: &[u8]) -> &mut Self {
        for b in bytes {
            self.0 ^= u64::from(*b);
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
        self
    }
    /// Length-prefixed, so concatenations cannot collide by shifting bytes.
    fn field(&mut self, text: &str) -> &mut Self {
        self.bytes(&(text.len() as u64).to_le_bytes())
            .bytes(text.as_bytes())
    }
    fn number(&mut self, n: u64) -> &mut Self {
        self.bytes(&n.to_le_bytes())
    }
    fn finish(&self) -> Fingerprint {
        Fingerprint(self.0)
    }
}

/// Fingerprints of the current inputs of every unit.
///
/// A unit fingerprint covers the package selection and result provenance,
/// the unit kind and specification, the feed composition for feeds, and for
/// every inlet port the fingerprint of the connected upstream unit and its
/// outlet port. Tags and stream names are not inputs. Units on a recycle or
/// downstream of an open inlet still receive fingerprints; units that cannot
/// be ordered because of a recycle receive none and therefore never match a
/// stored result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fingerprints {
    /// Package selection plus provenance, or `None` if the package cannot be built.
    pub package: Option<Fingerprint>,
    /// Per-unit fingerprints.
    pub units: BTreeMap<UnitId, Fingerprint>,
}

impl Fingerprints {
    /// Computes fingerprints for the given inputs.
    pub fn compute(fs: &Flowsheet, feeds: &FeedFractions, package: &PackageSelection) -> Self {
        Self::compute_with(fs, feeds, Self::package(package))
    }

    /// Fingerprint of a package selection together with the provenance of
    /// results it would produce, or `None` if the package cannot be built.
    /// Building the model is comparatively expensive; interactive callers
    /// cache this value per selection and use [`Self::compute_with`].
    pub fn package(package: &PackageSelection) -> Option<Fingerprint> {
        package.provenance().ok().map(|prov| {
            let text = serde_json::to_string(&(package, prov)).expect("serializable selection");
            Fnv::new().field("package").field(&text).finish()
        })
    }

    /// [`Self::compute`] with a precomputed [`Self::package`] fingerprint.
    pub fn compute_with(
        fs: &Flowsheet,
        feeds: &FeedFractions,
        package_fp: Option<Fingerprint>,
    ) -> Self {
        let mut units = BTreeMap::new();
        if let Some(pfp) = package_fp {
            let mut visiting = BTreeSet::new();
            for u in fs.units() {
                unit_fp(fs, feeds, pfp, u.id, &mut units, &mut visiting);
            }
        }
        Self {
            package: package_fp,
            units,
        }
    }

    /// Fingerprint of a stream: its source unit's fingerprint and outlet port.
    pub fn stream(&self, fs: &Flowsheet, id: StreamId) -> Option<Fingerprint> {
        let s = fs.stream(id)?;
        let source = self.units.get(&s.from.unit)?;
        Some(stream_fp(*source, s.from.port))
    }

    /// One digest of all unit fingerprints, used to recognise whether the
    /// whole input set changed since a solve attempt.
    pub fn combined(&self) -> Fingerprint {
        let mut h = Fnv::new();
        h.field("project");
        if let Some(p) = self.package {
            h.number(p.0);
        }
        for (id, fp) in &self.units {
            h.number(id.0).number(fp.0);
        }
        h.finish()
    }
}

fn stream_fp(source: Fingerprint, port: usize) -> Fingerprint {
    Fnv::new()
        .field("stream")
        .number(source.0)
        .number(port as u64)
        .finish()
}

fn unit_fp(
    fs: &Flowsheet,
    feeds: &FeedFractions,
    package: Fingerprint,
    id: UnitId,
    memo: &mut BTreeMap<UnitId, Fingerprint>,
    visiting: &mut BTreeSet<UnitId>,
) -> Option<Fingerprint> {
    if let Some(fp) = memo.get(&id) {
        return Some(*fp);
    }
    if !visiting.insert(id) {
        return None; // recycle
    }
    let unit = fs.unit(id)?;
    let mut h = Fnv::new();
    h.field("unit").number(package.0);
    h.field(&serde_json::to_string(&KindRecord::from(unit.kind)).expect("kind"));
    h.field(&serde_json::to_string(&SpecRecord::from(unit.spec)).expect("spec"));
    match feeds.get(&id) {
        Some(x) => h.field(&serde_json::to_string(x).expect("fractions")),
        None => h.field("no composition"),
    };
    let mut complete = true;
    for (i, port) in unit.kind.ports().iter().enumerate() {
        if port.direction != PortDirection::Inlet {
            continue;
        }
        h.number(i as u64);
        match fs.stream_at(PortRef::new(id, i)) {
            Some(s) => match unit_fp(fs, feeds, package, s.from.unit, memo, visiting) {
                Some(up) => {
                    h.number(stream_fp(up, s.from.port).0);
                }
                None => complete = false,
            },
            None => {
                h.field("open");
            }
        }
    }
    visiting.remove(&id);
    if !complete {
        return None;
    }
    let fp = h.finish();
    memo.insert(id, fp);
    Some(fp)
}
