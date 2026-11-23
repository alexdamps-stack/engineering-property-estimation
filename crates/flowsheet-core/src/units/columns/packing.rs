//! Packed column sections on the ideal theoretical-stage path.
//!
//! A packed section is represented in the equilibrium-stage solver by a
//! whole number of ideal theoretical stages, derived from the packed height
//! \(Z\) and the height equivalent to a theoretical plate (HETP):
//!
//! \[ N_\text{raw} = Z / \text{HETP}, \qquad N = \operatorname{round}(N_\text{raw}). \]
//!
//! The HETP concept and \(Z = N \cdot \text{HETP}\) follow Kister,
//! *Distillation Design* (McGraw-Hill, 1992), ch. 11, and Seader, Henley
//! and Roper, *Separation Process Principles*, 3rd ed. (Wiley, 2011),
//! §6.8. The HETP lumps mass-transfer efficiency into the height, so the
//! stages are ideal (Murphree efficiency 1); a stage efficiency is never
//! applied on top of an HETP.
//!
//! Two mutually exclusive bases ([`PackedSectionBasis`]):
//!
//! - **Packed height** (rating an existing bed): \(N\) follows from
//!   \(Z/\text{HETP}\) and the [`StageRounding`] rule. The solver models
//!   \(N\) stages, which correspond to an effective height \(N\cdot\text{HETP}\);
//!   the deviation from \(Z\) and the effective HETP \(Z/N\) are reported.
//! - **Theoretical stages** (design, as in RadFrac/DWSIM where the stage
//!   count is primary): \(N\) is given and the packed height
//!   \(Z = N\cdot\text{HETP}\) is derived exactly.
//!
//! Stage resolution: a quotient within [`INTEGER_SNAP_TOLERANCE`]
//! (relative) of a whole number is taken as that whole number before
//! rounding, so floating-point quotients such as `1.2 / 0.4 =
//! 2.9999999999999996` resolve to 3 under every rule. A resolved count
//! below 1 or above [`MAX_PACKED_STAGES`] fails explicitly; nothing is
//! clamped.
//!
//! HETP is a user input here. HETP correlations, mass-transfer and
//! hydraulic (pressure drop, flooding) correlations belong to the
//! rate-based packed model, which uses the same [`PackedSectionSpec`]
//! height but discretises it into segments instead of theoretical stages
//! (`docs/a17-columns-packing.md`).

use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional};

/// Largest number of theoretical stages one packed section may resolve to.
pub const MAX_PACKED_STAGES: u16 = 500;

/// Relative distance from a whole number below which `Z / HETP` is taken
/// as exactly that whole number before rounding.
pub const INTEGER_SNAP_TOLERANCE: f64 = 1e-9;

/// How a fractional `Z / HETP` becomes a whole number of stages.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StageRounding {
    /// Nearest whole number, halves rounded up. Smallest height deviation
    /// (at most half an HETP).
    #[default]
    Nearest,
    /// Round down: never credits the bed with more separation than its
    /// height provides (conservative when rating a given bed).
    Down,
    /// Round up: never models fewer stages than the height suggests.
    Up,
}

/// What determines the stage count of a packed section (mutually
/// exclusive).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PackedSectionBasis {
    /// Packed bed height in m (strictly positive); stages follow from
    /// `height / HETP` with the given rounding rule.
    PackedHeight {
        /// Packed height in m.
        height_m: Option<f64>,
        /// Rounding of the fractional stage count.
        rounding: StageRounding,
    },
    /// Number of theoretical stages, 1 to [`MAX_PACKED_STAGES`]; the packed
    /// height is `stages × HETP`.
    TheoreticalStages(u16),
}

impl Default for PackedSectionBasis {
    fn default() -> Self {
        PackedSectionBasis::PackedHeight {
            height_m: None,
            rounding: StageRounding::default(),
        }
    }
}

/// Specification of one packed column section on the equilibrium
/// (ideal theoretical-stage) path.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PackedSectionSpec {
    /// Height equivalent to a theoretical plate in m (strictly positive,
    /// required).
    pub hetp_m: Option<f64>,
    /// Height or stage basis.
    pub basis: PackedSectionBasis,
}

/// A packed section resolved to ideal theoretical stages, the input the
/// equilibrium-stage solver consumes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PackedStages {
    /// Number of ideal theoretical stages, 1 to [`MAX_PACKED_STAGES`].
    pub stages: u16,
    /// HETP in m as specified.
    pub hetp_m: f64,
    /// Packed height in m: the specified height, or `stages × HETP` on the
    /// stage basis.
    pub packed_height_m: f64,
    /// Unrounded `height / HETP` (equal to `stages` on the stage basis).
    pub raw_stages: f64,
}

impl PackedStages {
    /// Height the modelled stages represent, `stages × HETP`, in m.
    pub fn modelled_height_m(&self) -> f64 {
        f64::from(self.stages) * self.hetp_m
    }

    /// Modelled minus packed height in m: positive when rounding credits
    /// the bed with more stages than its height provides.
    pub fn height_deviation_m(&self) -> f64 {
        self.modelled_height_m() - self.packed_height_m
    }

    /// HETP implied by the whole stage count, `height / stages`, in m.
    pub fn effective_hetp_m(&self) -> f64 {
        self.packed_height_m / f64::from(self.stages)
    }

    /// Murphree vapour efficiency of every packed stage. Always 1: the
    /// HETP already accounts for non-ideal contacting.
    pub fn stage_efficiency(&self) -> f64 {
        1.0
    }

    /// Zero-based stage (counted from the top of the section) that holds
    /// the packing at `depth_m` below the top of the bed, for placing a feed
    /// or draw given by elevation. Stage `k` spans depths
    /// `[k, k + 1) × effective HETP`; the bottom of the bed belongs to the
    /// last stage. Fails outside `[0, packed height]`.
    pub fn stage_at_depth(&self, depth_m: f64) -> Result<u16, String> {
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.packed_height_m {
            return Err(format!(
                "depth {depth_m} m lies outside the packed bed (0 to {} m)",
                self.packed_height_m
            ));
        }
        let position = snap(depth_m / self.effective_hetp_m()).floor();
        // `position` is finite and within [0, stages]; the bottom face maps
        // to the last stage.
        Ok((position as u16).min(self.stages - 1))
    }
}

/// Snaps `x` to the nearest whole number if it is within
/// [`INTEGER_SNAP_TOLERANCE`] of it (relative, absolute near zero).
fn snap(x: f64) -> f64 {
    let whole = x.round();
    if (x - whole).abs() <= INTEGER_SNAP_TOLERANCE * whole.abs().max(1.0) {
        whole
    } else {
        x
    }
}

impl PackedSectionSpec {
    /// Missing and invalid fields, in the group-spec convention.
    pub fn issues(&self) -> Vec<SpecIssue> {
        let mut out = Vec::new();
        check_optional(&mut out, "HETP", self.hetp_m, Sign::Positive);
        match self.basis {
            PackedSectionBasis::PackedHeight { height_m, .. } => {
                check_optional(&mut out, "packed height", height_m, Sign::Positive)
            }
            PackedSectionBasis::TheoreticalStages(n) => {
                if !(1..=MAX_PACKED_STAGES).contains(&n) {
                    out.push(SpecIssue {
                        field: "theoretical stages",
                        problem: SpecProblem::CountOutOfRange,
                    });
                }
            }
        }
        out
    }

    /// Resolves the section to ideal theoretical stages. Fails on any
    /// specification issue, on a height shorter than the rounding rule
    /// allows (zero stages), and on more than [`MAX_PACKED_STAGES`].
    pub fn resolve(&self) -> Result<PackedStages, String> {
        if let Some(issue) = self.issues().first() {
            return Err(format!("packed section: {} {}", issue.field, issue.problem));
        }
        // Both unwraps are guarded by `issues`.
        let hetp_m = self.hetp_m.expect("HETP checked");
        match self.basis {
            PackedSectionBasis::TheoreticalStages(stages) => Ok(PackedStages {
                stages,
                hetp_m,
                packed_height_m: f64::from(stages) * hetp_m,
                raw_stages: f64::from(stages),
            }),
            PackedSectionBasis::PackedHeight { height_m, rounding } => {
                let height_m = height_m.expect("height checked");
                let raw_stages = height_m / hetp_m;
                if !raw_stages.is_finite() {
                    return Err(format!(
                        "packed section: height {height_m} m / HETP {hetp_m} m is not finite"
                    ));
                }
                let snapped = snap(raw_stages);
                let whole = match rounding {
                    StageRounding::Nearest => snapped.round(),
                    StageRounding::Down => snapped.floor(),
                    StageRounding::Up => snapped.ceil(),
                };
                if whole < 1.0 {
                    return Err(format!(
                        "packed section: height {height_m} m is {raw_stages:.4} HETP, \
                         which rounds to no theoretical stage"
                    ));
                }
                if whole > f64::from(MAX_PACKED_STAGES) {
                    return Err(format!(
                        "packed section: height {height_m} m is {raw_stages:.1} HETP, \
                         more than {MAX_PACKED_STAGES} theoretical stages"
                    ));
                }
                Ok(PackedStages {
                    stages: whole as u16,
                    hetp_m,
                    packed_height_m: height_m,
                    raw_stages,
                })
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn by_height(height: f64, hetp: f64, rounding: StageRounding) -> PackedSectionSpec {
        PackedSectionSpec {
            hetp_m: Some(hetp),
            basis: PackedSectionBasis::PackedHeight {
                height_m: Some(height),
                rounding,
            },
        }
    }

    fn stages(height: f64, hetp: f64, rounding: StageRounding) -> u16 {
        by_height(height, hetp, rounding).resolve().unwrap().stages
    }

    #[test]
    fn default_spec_reports_missing_hetp_and_height() {
        let issues = PackedSectionSpec::default().issues();
        assert_eq!(
            issues,
            [
                SpecIssue {
                    field: "HETP",
                    problem: SpecProblem::Missing
                },
                SpecIssue {
                    field: "packed height",
                    problem: SpecProblem::Missing
                },
            ]
        );
        let err = PackedSectionSpec::default().resolve().unwrap_err();
        assert_eq!(err, "packed section: HETP is missing");
    }

    #[test]
    fn nonpositive_and_nonfinite_inputs_are_rejected() {
        for bad in [0.0, -0.5, f64::NAN, f64::INFINITY] {
            let hetp = by_height(5.0, bad, StageRounding::Nearest).issues();
            assert_eq!(hetp.len(), 1, "HETP {bad}");
            assert_eq!(hetp[0].field, "HETP");
            let height = by_height(bad, 0.5, StageRounding::Nearest).issues();
            assert_eq!(height.len(), 1, "height {bad}");
            assert_eq!(height[0].field, "packed height");
        }
        assert_eq!(
            by_height(-1.0, 0.5, StageRounding::Nearest).issues()[0].problem,
            SpecProblem::NotPositive
        );
        assert_eq!(
            by_height(f64::NAN, 0.5, StageRounding::Nearest).issues()[0].problem,
            SpecProblem::NotFinite
        );
    }

    #[test]
    fn representative_bed_resolves_with_each_rounding_rule() {
        // 6.2 m of packing with HETP 0.45 m: 13.78 theoretical stages.
        let s = by_height(6.2, 0.45, StageRounding::Nearest)
            .resolve()
            .unwrap();
        assert!((s.raw_stages - 6.2 / 0.45).abs() < 1e-12);
        assert_eq!(s.stages, 14);
        assert!((s.modelled_height_m() - 6.3).abs() < 1e-12);
        assert!((s.height_deviation_m() - 0.1).abs() < 1e-12);
        assert!((s.effective_hetp_m() - 6.2 / 14.0).abs() < 1e-15);
        assert_eq!(s.stage_efficiency(), 1.0);
        assert_eq!(stages(6.2, 0.45, StageRounding::Down), 13);
        assert_eq!(stages(6.2, 0.45, StageRounding::Up), 14);
        // Halves round up under Nearest.
        assert_eq!(stages(2.5, 1.0, StageRounding::Nearest), 3);
        assert_eq!(stages(2.49, 1.0, StageRounding::Nearest), 2);
    }

    #[test]
    fn exact_multiples_are_not_perturbed_by_floating_point() {
        // 1.2 / 0.4 = 2.9999999999999996 and 0.7 / 0.1 = 6.999999999999999.
        assert!(std::hint::black_box(1.2_f64) / 0.4 < 3.0);
        for rounding in [
            StageRounding::Nearest,
            StageRounding::Down,
            StageRounding::Up,
        ] {
            assert_eq!(stages(1.2, 0.4, rounding), 3);
            assert_eq!(stages(0.7, 0.1, rounding), 7);
            assert_eq!(stages(0.3, 0.1, rounding), 3);
            assert_eq!(stages(10.0, 0.5, rounding), 20);
        }
        // Genuinely fractional values are not snapped.
        assert_eq!(stages(3.0 + 1e-6, 1.0, StageRounding::Up), 4);
        assert_eq!(stages(3.0 - 1e-6, 1.0, StageRounding::Down), 2);
    }

    #[test]
    fn short_beds_fail_instead_of_producing_zero_stages() {
        let err = by_height(0.3, 0.5, StageRounding::Down)
            .resolve()
            .unwrap_err();
        assert!(err.contains("rounds to no theoretical stage"), "{err}");
        assert!(
            by_height(0.2, 0.5, StageRounding::Nearest)
                .resolve()
                .is_err()
        );
        assert_eq!(stages(0.25, 0.5, StageRounding::Nearest), 1);
        assert_eq!(stages(0.01, 0.5, StageRounding::Up), 1);
        // One HETP exactly is one stage under every rule.
        assert_eq!(stages(0.5, 0.5, StageRounding::Down), 1);
    }

    #[test]
    fn stage_count_is_capped_without_clamping() {
        let max = f64::from(MAX_PACKED_STAGES);
        assert_eq!(
            stages(max * 0.2, 0.2, StageRounding::Nearest),
            MAX_PACKED_STAGES
        );
        let err = by_height((max + 1.0) * 0.2, 0.2, StageRounding::Nearest)
            .resolve()
            .unwrap_err();
        assert!(err.contains("more than 500"), "{err}");
        // Overflowing quotient: tiny HETP with a large height.
        assert!(
            by_height(1e300, 1e-300, StageRounding::Down)
                .resolve()
                .is_err()
        );
    }

    #[test]
    fn stage_basis_derives_height_exactly() {
        let spec = PackedSectionSpec {
            hetp_m: Some(0.6),
            basis: PackedSectionBasis::TheoreticalStages(8),
        };
        let s = spec.resolve().unwrap();
        assert_eq!(s.stages, 8);
        assert!((s.packed_height_m - 4.8).abs() < 1e-12);
        assert_eq!(s.height_deviation_m(), 0.0);
        assert_eq!(s.effective_hetp_m(), 0.6);
        for n in [0, MAX_PACKED_STAGES + 1] {
            let bad = PackedSectionSpec {
                hetp_m: Some(0.6),
                basis: PackedSectionBasis::TheoreticalStages(n),
            };
            assert_eq!(
                bad.issues(),
                [SpecIssue {
                    field: "theoretical stages",
                    problem: SpecProblem::CountOutOfRange
                }]
            );
            assert!(bad.resolve().is_err());
        }
    }

    #[test]
    fn depth_maps_to_stages_from_the_top() {
        // 4 m bed, HETP 1 m: four stages of 1 m each.
        let s = by_height(4.0, 1.0, StageRounding::Nearest)
            .resolve()
            .unwrap();
        assert_eq!(s.stage_at_depth(0.0), Ok(0));
        assert_eq!(s.stage_at_depth(0.999), Ok(0));
        assert_eq!(s.stage_at_depth(1.0), Ok(1));
        assert_eq!(s.stage_at_depth(2.5), Ok(2));
        assert_eq!(s.stage_at_depth(4.0), Ok(3));
        assert!(s.stage_at_depth(-0.1).is_err());
        assert!(s.stage_at_depth(4.1).is_err());
        assert!(s.stage_at_depth(f64::NAN).is_err());
        // Rounded bed: 6.2 m as 14 stages of 6.2/14 m; depth uses the
        // effective HETP so the whole bed maps onto the modelled stages.
        let r = by_height(6.2, 0.45, StageRounding::Nearest)
            .resolve()
            .unwrap();
        assert_eq!(r.stage_at_depth(6.2), Ok(13));
        assert_eq!(r.stage_at_depth(3.1), Ok(7));
        // Stage boundary computed in floating point.
        let t = by_height(1.2, 0.4, StageRounding::Nearest)
            .resolve()
            .unwrap();
        assert_eq!(t.stage_at_depth(0.8), Ok(2));
    }
}
