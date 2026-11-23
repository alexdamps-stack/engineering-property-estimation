//! Reactor chemistry editors: stoichiometric reactions of a conversion
//! reactor and mass yields of a yield reactor (#191, #192).
//!
//! Coefficients and yields are edited per package component; a zero entry is
//! not stored. Components named in the specification but missing from the
//! package (for example after the package changed) stay listed so they can be
//! removed; the engine reports them before solving.

use eframe::egui::{self, ComboBox, DragValue, RichText, Ui};
use flowsheet_core::{
    Arrhenius, ConcentrationBasis, Correlation, EquilibriumConstant, EquilibriumReaction,
    KineticReaction, Kinetics, LhhwTerm, RateBasis, RateLaw, RateTerm, ReactingPhase, Reaction,
    ReactionExtent, UnitSpec,
};

use super::spec::{choice, opt_field};
use crate::symbols::Theme;
use crate::units::{DisplayUnits, Quantity};

/// Chemistry section of a reactor; does nothing for other specifications.
pub(super) fn chemistry_form(
    ui: &mut Ui,
    spec: &mut UnitSpec,
    components: &[String],
    u: DisplayUnits,
) {
    match spec {
        UnitSpec::ConversionReactor { reactions, .. } => {
            reactions_form(ui, reactions, components, u)
        }
        UnitSpec::YieldReactor {
            mass_yields,
            inerts,
            ..
        } => yields_form(ui, mass_yields, inerts, components),
        UnitSpec::EquilibriumReactor { reactions, .. } => {
            equilibrium_form(ui, reactions, components)
        }
        UnitSpec::GibbsReactor { inerts, .. } => inerts_form(ui, inerts, components),
        UnitSpec::Cstr { kinetics, .. }
        | UnitSpec::PlugFlow { kinetics, .. }
        | UnitSpec::Batch { kinetics, .. } => kinetics_form(ui, kinetics, components),
        _ => {}
    }
}

/// Package components followed by any other identifiers the entries name.
fn listed<'a>(
    components: &'a [String],
    named: impl Iterator<Item = &'a String>,
) -> Vec<&'a String> {
    let mut all: Vec<&String> = components.iter().collect();
    for id in named {
        if !all.contains(&id) {
            all.push(id);
        }
    }
    all
}

/// Sets `id` to `value` in an ordered `(id, value)` list; zero removes it.
fn set_entry(entries: &mut Vec<(String, f64)>, id: &str, value: f64) {
    match entries.iter().position(|(c, _)| c == id) {
        Some(i) if value == 0.0 => {
            entries.remove(i);
        }
        Some(i) => entries[i].1 = value,
        None if value != 0.0 => entries.push((id.to_owned(), value)),
        None => {}
    }
}

/// One grid row per package component (plus any unknown names) editing a
/// sparse `(component, value)` list; zero removes the entry.
fn entry_rows(
    ui: &mut Ui,
    entries: &mut Vec<(String, f64)>,
    components: &[String],
    theme: &Theme,
    hover: &str,
) {
    let named: Vec<String> = entries.iter().map(|(c, _)| c.clone()).collect();
    for id in listed(components, named.iter()) {
        if components.contains(id) {
            ui.label(id.as_str());
        } else {
            ui.colored_label(theme.warning, format!("{id} (not in package)"));
        }
        let mut v = entries
            .iter()
            .find(|(c, _)| c == id)
            .map_or(0.0, |(_, v)| *v);
        if ui
            .add(DragValue::new(&mut v).speed(0.1).max_decimals(6))
            .on_hover_text(hover)
            .changed()
        {
            set_entry(entries, id, v);
        }
        ui.end_row();
    }
}

/// `2 methane + oxygen -> ...` with unit coefficients omitted.
pub(crate) fn equation(coefficients: &[(String, f64)]) -> String {
    let side = |sign: f64| {
        coefficients
            .iter()
            .filter(|(_, v)| v.is_finite() && *v * sign > 0.0)
            .map(|(id, v)| {
                let v = v.abs();
                if v == 1.0 {
                    id.clone()
                } else {
                    format!("{} {id}", trim(v))
                }
            })
            .collect::<Vec<_>>()
            .join(" + ")
    };
    format!("{} → {}", side(-1.0), side(1.0))
}

fn trim(v: f64) -> String {
    let s = format!("{v:.6}");
    s.trim_end_matches('0').trim_end_matches('.').to_owned()
}

fn reactions_form(
    ui: &mut Ui,
    reactions: &mut Vec<Reaction>,
    components: &[String],
    u: DisplayUnits,
) {
    ui.add_space(8.0);
    ui.strong("Reactions");
    ui.label(
        RichText::new(
            "Applied in order. Coefficients: negative for reactants, positive for products; every element must balance.",
        )
        .small()
        .weak(),
    );
    let theme = Theme::from_visuals(ui.visuals());
    let mut remove = None;
    for (k, r) in reactions.iter_mut().enumerate() {
        ui.add_space(4.0);
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(RichText::new(format!("Reaction {}", k + 1)).strong());
                if ui.small_button("Remove").clicked() {
                    remove = Some(k);
                }
            });
            ui.label(RichText::new(equation(&r.coefficients)).monospace());
            egui::Grid::new(("coefficients", k))
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    entry_rows(
                        ui,
                        &mut r.coefficients,
                        components,
                        &theme,
                        "Stoichiometric coefficient; 0 = not taking part",
                    );
                    extent_form(ui, &mut r.extent, &r.coefficients, u);
                });
        });
    }
    if let Some(k) = remove {
        reactions.remove(k);
    }
    if ui.button("Add reaction").clicked() {
        reactions.push(Reaction {
            coefficients: Vec::new(),
            extent: ReactionExtent::Conversion {
                key: None,
                fraction: None,
            },
        });
    }
}

fn extent_form(
    ui: &mut Ui,
    extent: &mut ReactionExtent,
    coefficients: &[(String, f64)],
    u: DisplayUnits,
) {
    let current = matches!(extent, ReactionExtent::MolarExtent(_)) as usize;
    match choice(ui, "Extent", current, &["Conversion", "Molar extent"]) {
        Some(0) => {
            *extent = ReactionExtent::Conversion {
                key: None,
                fraction: None,
            }
        }
        Some(_) => *extent = ReactionExtent::MolarExtent(None),
        None => {}
    }
    match extent {
        ReactionExtent::Conversion { key, fraction } => {
            let label = ui.label("Key reactant");
            let reactants: Vec<&String> = coefficients
                .iter()
                .filter(|(_, v)| *v < 0.0)
                .map(|(c, _)| c)
                .collect();
            // The id is scoped by the enclosing per-reaction grid.
            ComboBox::from_id_salt("key")
                .selected_text(key.as_deref().unwrap_or("choose"))
                .show_ui(ui, |ui| {
                    for id in reactants {
                        ui.selectable_value(key, Some(id.clone()), id.as_str());
                    }
                })
                .response
                .labelled_by(label.id);
            ui.end_row();
            opt_field(ui, "Conversion", fraction, Quantity::Fraction, u);
        }
        ReactionExtent::MolarExtent(e) => opt_field(ui, "Extent", e, Quantity::MolarFlow, u),
    }
}

fn yields_form(
    ui: &mut Ui,
    yields: &mut Vec<(String, f64)>,
    inerts: &mut Vec<String>,
    components: &[String],
) {
    ui.add_space(8.0);
    ui.strong("Mass yields");
    ui.label(
        RichText::new(
            "Relative weights of the outlet mass from the non-inert feed; normalized so mass is conserved. Inert components pass through.",
        )
        .small()
        .weak(),
    );
    let theme = Theme::from_visuals(ui.visuals());
    let named: Vec<String> = yields
        .iter()
        .map(|(c, _)| c.clone())
        .chain(inerts.iter().cloned())
        .collect();
    let total: f64 = yields.iter().map(|(_, y)| y.max(0.0)).sum();
    egui::Grid::new("yields")
        .num_columns(4)
        .spacing([8.0, 4.0])
        .show(ui, |ui| {
            ui.label("");
            ui.label(RichText::new("weight").weak());
            ui.label(RichText::new("share").weak());
            ui.label(RichText::new("inert").weak());
            ui.end_row();
            for id in listed(components, named.iter()) {
                let known = components.contains(id);
                if known {
                    ui.label(id.as_str());
                } else {
                    ui.colored_label(theme.warning, format!("{id} (not in package)"));
                }
                let mut inert = inerts.contains(id);
                let mut y = yields
                    .iter()
                    .find(|(c, _)| c == id)
                    .map_or(0.0, |(_, y)| *y);
                let edited = ui
                    .add_enabled(
                        !inert,
                        DragValue::new(&mut y)
                            .speed(0.05)
                            .range(0.0..=f64::MAX)
                            .clamp_existing_to_range(false)
                            .max_decimals(6),
                    )
                    .changed();
                if edited {
                    set_entry(yields, id, y);
                }
                ui.label(if total > 0.0 && !inert {
                    format!("{:.1} %", 100.0 * y.max(0.0) / total)
                } else {
                    String::new()
                });
                if ui.checkbox(&mut inert, "").changed() {
                    inerts.retain(|c| c != id);
                    if inert {
                        inerts.push(id.clone());
                        set_entry(yields, id, 0.0);
                    }
                }
                ui.end_row();
            }
        });
}

/// Frame with a numbered title and a remove button; returns whether the
/// item should be removed.
fn item_frame(ui: &mut Ui, title: &str, body: impl FnOnce(&mut Ui)) -> bool {
    let mut remove = false;
    ui.add_space(4.0);
    egui::Frame::group(ui.style()).show(ui, |ui| {
        ui.horizontal(|ui| {
            ui.label(RichText::new(title).strong());
            remove = ui.small_button("Remove").clicked();
        });
        body(ui);
    });
    remove
}

/// A drag value shown and typed in scientific notation, for constants that
/// span many orders of magnitude.
fn scientific(v: &mut f64) -> DragValue<'_> {
    let speed = (v.abs() * 0.01).max(1e-12);
    DragValue::new(v)
        .speed(speed)
        .custom_formatter(|x, _| format!("{x:.6e}"))
        .custom_parser(|s| s.trim().parse().ok())
}

/// `a + b/T + c ln T + d T` coefficients in one row.
fn correlation_row(ui: &mut Ui, label: &str, c: &mut Correlation) {
    ui.label(label)
        .on_hover_text("ln K = a + b/T + c ln T + d T, T in K");
    ui.horizontal(|ui| {
        for (name, v) in [
            ("a", &mut c.a),
            ("b", &mut c.b),
            ("c", &mut c.c),
            ("d", &mut c.d),
        ] {
            ui.add(scientific(v).prefix(format!("{name} ")));
        }
    });
    ui.end_row();
}

fn equilibrium_form(ui: &mut Ui, reactions: &mut Vec<EquilibriumReaction>, components: &[String]) {
    ui.add_space(8.0);
    ui.strong("Equilibrium reactions");
    ui.label(
        RichText::new(
            "Independent reactions; every element must balance. K on the fugacity basis, standard state ideal gas at 1 bar.",
        )
        .small()
        .weak(),
    );
    let theme = Theme::from_visuals(ui.visuals());
    let mut remove = None;
    for (k, r) in reactions.iter_mut().enumerate() {
        if item_frame(ui, &format!("Reaction {}", k + 1), |ui| {
            ui.label(RichText::new(equation(&r.coefficients)).monospace());
            egui::Grid::new(("equilibrium", k))
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    entry_rows(
                        ui,
                        &mut r.coefficients,
                        components,
                        &theme,
                        "Stoichiometric coefficient; 0 = not taking part",
                    );
                    let current =
                        matches!(r.constant, EquilibriumConstant::Correlation(_)) as usize;
                    match choice(ui, "K from", current, &["Gibbs energies", "Correlation"]) {
                        Some(0) => r.constant = EquilibriumConstant::GibbsEnergy,
                        Some(_) => {
                            r.constant = EquilibriumConstant::Correlation(Correlation::default())
                        }
                        None => {}
                    }
                    if let EquilibriumConstant::Correlation(c) = &mut r.constant {
                        correlation_row(ui, "ln K", c);
                    }
                });
        }) {
            remove = Some(k);
        }
    }
    if let Some(k) = remove {
        reactions.remove(k);
    }
    if ui.button("Add reaction").clicked() {
        reactions.push(EquilibriumReaction {
            coefficients: Vec::new(),
            constant: EquilibriumConstant::GibbsEnergy,
        });
    }
}

fn inerts_form(ui: &mut Ui, inerts: &mut Vec<String>, components: &[String]) {
    ui.add_space(8.0);
    ui.strong("Inert components");
    ui.label(
        RichText::new("All other components may form or react; elements are conserved.")
            .small()
            .weak(),
    );
    let named = inerts.clone();
    for id in listed(components, named.iter()) {
        let mut inert = inerts.contains(id);
        let label = if components.contains(id) {
            id.clone()
        } else {
            format!("{id} (not in package)")
        };
        if ui.checkbox(&mut inert, label).changed() {
            inerts.retain(|c| c != id);
            if inert {
                inerts.push(id.clone());
            }
        }
    }
}

/// Arrhenius parameters as grid rows.
fn arrhenius_rows(ui: &mut Ui, k: &mut Arrhenius) {
    let k0 = ui.label("k0").on_hover_text(
        "Pre-exponential factor in the units implied by the basis, the orders and the rate basis",
    );
    ui.horizontal(|ui| match &mut k.pre_exponential {
        Some(v) => {
            ui.add(scientific(v)).labelled_by(k0.id);
            if ui.small_button("×").on_hover_text("Clear").clicked() {
                k.pre_exponential = None;
            }
        }
        None => {
            let warning = Theme::from_visuals(ui.visuals()).warning;
            ui.label(RichText::new("required").italics().color(warning));
            if ui.small_button("Enter value").clicked() {
                k.pre_exponential = Some(1.0);
            }
        }
    });
    ui.end_row();
    let e = ui.label("E");
    ui.add(
        DragValue::new(&mut k.activation_energy_j_mol)
            .speed(100.0)
            .suffix(" J/mol"),
    )
    .labelled_by(e.id);
    ui.end_row();
    let n = ui.label("n").on_hover_text("Temperature exponent");
    ui.add(DragValue::new(&mut k.temperature_exponent).speed(0.01))
        .labelled_by(n.id);
    ui.end_row();
    let mut referenced = k.reference_temperature_k.is_some();
    let t_ref = ui.label("T ref");
    ui.horizontal(|ui| {
        if ui
            .checkbox(&mut referenced, "")
            .labelled_by(t_ref.id)
            .on_hover_text(
                "k0 at a reference temperature: k = k0 (T/Tref)^n exp(-E/R (1/T - 1/Tref))",
            )
            .changed()
        {
            k.reference_temperature_k = referenced.then_some(298.15);
        }
        if let Some(t) = &mut k.reference_temperature_k {
            ui.add(DragValue::new(t).speed(1.0).suffix(" K"))
                .labelled_by(t_ref.id);
        } else {
            ui.label(RichText::new("k = k0 T^n exp(-E/RT)").small().weak());
        }
    });
    ui.end_row();
}

fn lhhw_term(ui: &mut Ui, label: &str, term: &mut LhhwTerm, components: &[String], theme: &Theme) {
    correlation_row(ui, label, &mut term.constant);
    entry_rows(
        ui,
        &mut term.orders,
        components,
        theme,
        "Order; 0 = not in this term",
    );
}

const BASES: [(ConcentrationBasis, &str); 6] = [
    (
        ConcentrationBasis::MolarConcentration,
        "Molar concentration (mol/m³)",
    ),
    (
        ConcentrationBasis::MassConcentration,
        "Mass concentration (kg/m³)",
    ),
    (ConcentrationBasis::MoleFraction, "Mole fraction"),
    (ConcentrationBasis::MassFraction, "Mass fraction"),
    (ConcentrationBasis::PartialPressure, "Partial pressure (Pa)"),
    (ConcentrationBasis::Fugacity, "Fugacity (Pa)"),
];

fn kinetics_form(ui: &mut Ui, kinetics: &mut Kinetics, components: &[String]) {
    ui.add_space(8.0);
    ui.strong("Kinetic reactions");
    ui.label(
        RichText::new(
            "Rates add. Rate per reacting-phase volume (mol/(m³ s)) or per catalyst mass (mol/(kg s)); every element must balance.",
        )
        .small()
        .weak(),
    );
    let theme = Theme::from_visuals(ui.visuals());
    egui::Grid::new("catalyst").num_columns(2).show(ui, |ui| {
        ui.label("Catalyst bulk density");
        ui.horizontal(|ui| match &mut kinetics.catalyst_bulk_density_kg_m3 {
            Some(v) => {
                ui.add(
                    DragValue::new(v)
                        .speed(1.0)
                        .range(0.0..=f64::MAX)
                        .suffix(" kg/m³"),
                );
                if ui.small_button("×").on_hover_text("Clear").clicked() {
                    kinetics.catalyst_bulk_density_kg_m3 = None;
                }
            }
            None => {
                ui.label(RichText::new("none").weak());
                if ui.small_button("Enter value").clicked() {
                    kinetics.catalyst_bulk_density_kg_m3 = Some(1000.0);
                }
            }
        });
        ui.end_row();
    });
    let mut remove = None;
    for (k, r) in kinetics.reactions.iter_mut().enumerate() {
        if item_frame(ui, &format!("Reaction {}", k + 1), |ui| {
            ui.label(RichText::new(equation(&r.coefficients)).monospace());
            egui::Grid::new(("kinetic", k))
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| kinetic_rows(ui, r, components, &theme));
        }) {
            remove = Some(k);
        }
    }
    if let Some(k) = remove {
        kinetics.reactions.remove(k);
    }
    if ui.button("Add reaction").clicked() {
        kinetics.reactions.push(KineticReaction {
            coefficients: Vec::new(),
            phase: ReactingPhase::Vapor,
            basis: ConcentrationBasis::MolarConcentration,
            rate_basis: RateBasis::ReactorVolume,
            rate: RateLaw::PowerLaw {
                forward: RateTerm::default(),
                reverse: None,
            },
        });
    }
}

fn kinetic_rows(ui: &mut Ui, r: &mut KineticReaction, components: &[String], theme: &Theme) {
    entry_rows(
        ui,
        &mut r.coefficients,
        components,
        theme,
        "Stoichiometric coefficient; 0 = not taking part",
    );
    let phase = matches!(r.phase, ReactingPhase::Liquid) as usize;
    match choice(ui, "Phase", phase, &["Vapor", "Liquid"]) {
        Some(0) => r.phase = ReactingPhase::Vapor,
        Some(_) => r.phase = ReactingPhase::Liquid,
        None => {}
    }
    let label = ui.label("Basis");
    ComboBox::from_id_salt("basis")
        .selected_text(
            BASES
                .iter()
                .find(|(b, _)| *b == r.basis)
                .map_or("", |(_, l)| l),
        )
        .show_ui(ui, |ui| {
            for (b, l) in BASES {
                ui.selectable_value(&mut r.basis, b, l);
            }
        })
        .response
        .labelled_by(label.id);
    ui.end_row();
    let per = matches!(r.rate_basis, RateBasis::CatalystMass) as usize;
    match choice(ui, "Rate per", per, &["Reactor volume", "Catalyst mass"]) {
        Some(0) => r.rate_basis = RateBasis::ReactorVolume,
        Some(_) => r.rate_basis = RateBasis::CatalystMass,
        None => {}
    }
    let law = matches!(r.rate, RateLaw::Lhhw { .. }) as usize;
    match choice(ui, "Rate law", law, &["Power law", "LHHW"]) {
        Some(0) => {
            r.rate = RateLaw::PowerLaw {
                forward: RateTerm::default(),
                reverse: None,
            }
        }
        Some(_) => {
            r.rate = RateLaw::Lhhw {
                kinetic: Arrhenius::default(),
                driving_forward: LhhwTerm::default(),
                driving_reverse: None,
                adsorption: Vec::new(),
                exponent: 1.0,
            }
        }
        None => {}
    }
    match &mut r.rate {
        RateLaw::PowerLaw { forward, reverse } => {
            ui.label(RichText::new("Forward").strong());
            ui.end_row();
            arrhenius_rows(ui, &mut forward.constant);
            entry_rows(
                ui,
                &mut forward.orders,
                components,
                theme,
                "Order; 0 = not in the rate",
            );
            let mut reversible = reverse.is_some();
            ui.label("");
            if ui.checkbox(&mut reversible, "Reversible").changed() {
                *reverse = reversible.then(RateTerm::default);
            }
            ui.end_row();
            if let Some(rev) = reverse {
                ui.label(RichText::new("Reverse").strong());
                ui.end_row();
                arrhenius_rows(ui, &mut rev.constant);
                entry_rows(
                    ui,
                    &mut rev.orders,
                    components,
                    theme,
                    "Order; 0 = not in the rate",
                );
            }
        }
        RateLaw::Lhhw {
            kinetic,
            driving_forward,
            driving_reverse,
            adsorption,
            exponent,
        } => {
            ui.label(RichText::new("Kinetic factor").strong());
            ui.end_row();
            arrhenius_rows(ui, kinetic);
            ui.label(RichText::new("Driving force").strong());
            ui.end_row();
            lhhw_term(ui, "Forward K", driving_forward, components, theme);
            let mut reversible = driving_reverse.is_some();
            ui.label("");
            if ui.checkbox(&mut reversible, "Reverse term").changed() {
                *driving_reverse = reversible.then(LhhwTerm::default);
            }
            ui.end_row();
            if let Some(rev) = driving_reverse {
                lhhw_term(ui, "Reverse K", rev, components, theme);
            }
            ui.label(RichText::new("Adsorption").strong());
            ui.end_row();
            let mut remove = None;
            for (i, term) in adsorption.iter_mut().enumerate() {
                ui.label(format!("Term {}", i + 1));
                if ui.small_button("Remove term").clicked() {
                    remove = Some(i);
                }
                ui.end_row();
                lhhw_term(ui, "K", term, components, theme);
            }
            if let Some(i) = remove {
                adsorption.remove(i);
            }
            ui.label("");
            if ui.small_button("Add adsorption term").clicked() {
                adsorption.push(LhhwTerm::default());
            }
            ui.end_row();
            ui.label("Exponent m");
            ui.add(DragValue::new(exponent).speed(0.1).range(0.0..=f64::MAX));
            ui.end_row();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equations_omit_unit_coefficients() {
        let r = Reaction {
            coefficients: vec![
                ("methane".into(), -1.0),
                ("oxygen".into(), -2.0),
                ("carbon-dioxide".into(), 1.0),
                ("water".into(), 2.0),
            ],
            extent: ReactionExtent::MolarExtent(None),
        };
        assert_eq!(
            equation(&r.coefficients),
            "methane + 2 oxygen → carbon-dioxide + 2 water"
        );
        let half = Reaction {
            coefficients: vec![("ethane".into(), -1.0), ("oxygen".into(), -3.5)],
            ..r
        };
        assert_eq!(equation(&half.coefficients), "ethane + 3.5 oxygen → ");
    }

    #[test]
    fn zero_entries_are_not_stored() {
        let mut e = vec![("a".to_string(), 1.0)];
        set_entry(&mut e, "b", 0.0);
        set_entry(&mut e, "c", 2.0);
        set_entry(&mut e, "a", 0.0);
        assert_eq!(e, [("c".to_string(), 2.0)]);
    }
}
