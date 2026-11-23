//! Technical-drawing look: CAD colour scheme and process-flow-diagram symbols.
//!
//! Equipment symbols are the traced vector geometry of the ISO 10628-2:2012
//! symbol sheets by Steven Baltakatei Sandoval (Baltakatei), 2020, updated
//! 2023 (<https://commons.wikimedia.org/wiki/File:ISO_10628-2_2012_Symbols.pdf>),
//! licensed CC BY-SA 4.0 (<https://creativecommons.org/licenses/by-sa/4.0/>).
//! The standalone SVG files extracted in `docs/design-reference/iso-10628/
//! symbols/` are embedded unchanged and parsed with `usvg` on first use;
//! their `<desc>` keeps the attribution and license. Adaptations made here:
//! the pipe lines drawn beside each symbol on the sheet are left out (the
//! canvas streams take their place), the geometry is scaled uniformly into
//! the unit's block, curves are flattened to polylines, the line work is
//! recoloured for the theme and state, and the body is masked with the
//! background colour.
//!
//! [`TRACED`] lists the asset chosen for each unit kind. Connection ports
//! are not part of the drawings: they are defined in [`crate::view`], and a
//! test checks that each one lies on the traced line work. The splitter and
//! the feed/product flags have no fitting reference symbol and remain
//! project drawings, as does the heating/cooling medium arrow on the heater
//! and cooler. Arrowheads are closed and filled with a 15° half angle as in
//! ISO 128. The dark scheme mirrors a CAD model space; the light scheme a
//! plotted sheet.

use std::sync::OnceLock;

use eframe::egui::{Color32, Painter, Pos2, Rect, Shape, Stroke, Vec2, Visuals, pos2, vec2};
use flowsheet_core::UnitKind;
use usvg::tiny_skia_path::{PathSegment, Point};

/// Colours used on the canvas, derived from the active light/dark visuals.
#[derive(Debug, Clone, Copy)]
pub struct Theme {
    /// Drawing background (model space or paper).
    pub background: Color32,
    /// Minor grid lines.
    pub grid_minor: Color32,
    /// Major grid lines (every fifth minor line).
    pub grid_major: Color32,
    /// Line work and lettering.
    pub line: Color32,
    /// De-emphasised line work (underspecified units, hints).
    pub muted: Color32,
    /// Selection highlight and grips.
    pub select: Color32,
    /// Snap marker on free nozzles.
    pub snap: Color32,
    /// Warnings (amber), on the canvas and in the surrounding panels.
    pub warning: Color32,
    /// Failures (red).
    pub error: Color32,
    /// Current, solved results (green check badges).
    pub ok: Color32,
}

impl Theme {
    /// Theme for the given egui visuals.
    pub fn from_visuals(visuals: &Visuals) -> Self {
        if visuals.dark_mode {
            Self {
                background: Color32::from_rgb(33, 40, 48),
                grid_minor: Color32::from_rgb(40, 48, 57),
                grid_major: Color32::from_rgb(54, 64, 76),
                line: Color32::from_rgb(240, 240, 240),
                muted: Color32::from_rgb(128, 138, 150),
                select: Color32::from_rgb(30, 144, 255),
                snap: Color32::from_rgb(255, 200, 0),
                warning: Color32::from_rgb(210, 153, 34),
                error: Color32::from_rgb(248, 81, 73),
                ok: Color32::from_rgb(63, 185, 80),
            }
        } else {
            Self {
                background: Color32::WHITE,
                grid_minor: Color32::from_rgb(240, 242, 245),
                grid_major: Color32::from_rgb(222, 226, 232),
                line: Color32::BLACK,
                muted: Color32::from_rgb(140, 140, 140),
                select: Color32::from_rgb(0, 102, 214),
                snap: Color32::from_rgb(214, 130, 0),
                warning: Color32::from_rgb(154, 103, 0),
                error: Color32::from_rgb(207, 34, 46),
                ok: Color32::from_rgb(26, 127, 55),
            }
        }
    }
}

/// Closed, filled arrowhead with a 15° half angle (ISO 128), tip at `tip`
/// pointing along `dir`.
pub fn arrowhead(painter: &Painter, tip: Pos2, dir: Vec2, length: f32, color: Color32) {
    let d = dir.normalized();
    let n = d.rot90() * length * (15.0_f32).to_radians().tan();
    let base = tip - d * length;
    painter.add(Shape::convex_polygon(
        vec![tip, base + n, base - n],
        color,
        Stroke::NONE,
    ));
}

/// A unit symbol drawn from a reference SVG in
/// `docs/design-reference/iso-10628/symbols/`.
pub struct TracedSymbol {
    /// Unit kind drawn with this symbol.
    pub kind: UnitKind,
    /// Asset ID in `docs/design-reference/iso-10628/manifest.json`.
    pub asset: &'static str,
    /// The SVG file, embedded unchanged.
    svg: &'static str,
    /// IDs of the pipe lines drawn beside the symbol on the reference sheet.
    /// They show where lines connect and are left out on the canvas, where
    /// the streams take their place.
    pipes: &'static [&'static str],
    /// IDs of nozzle stubs that are drawn but are not part of the body, so
    /// the background mask does not cover the gap between stub and body.
    stubs: &'static [&'static str],
    /// Whether the kind's port coordinates in [`crate::view`] were reviewed
    /// against the traced line work (a test then checks every port lies on
    /// it). A17 kinds start unreviewed: their nozzles sit on the block edge
    /// and the owning group reviews them with the calculation.
    pub reviewed_ports: bool,
}

macro_rules! traced {
    ($kind:ident, $asset:literal, pipes: $pipes:expr, stubs: $stubs:expr $(,)?) => {
        TracedSymbol {
            kind: UnitKind::$kind,
            asset: $asset,
            svg: include_str!(concat!(
                "../../../docs/design-reference/iso-10628/symbols/",
                $asset,
                ".svg"
            )),
            pipes: &$pipes,
            stubs: &$stubs,
            reviewed_ports: true,
        }
    };
    // A17 kinds: the whole reference drawing, nozzle stubs included.
    ($kind:ident, $asset:literal $(,)?) => {
        TracedSymbol {
            kind: UnitKind::$kind,
            asset: $asset,
            svg: include_str!(concat!(
                "../../../docs/design-reference/iso-10628/symbols/",
                $asset,
                ".svg"
            )),
            pipes: &[],
            stubs: &[],
            reviewed_ports: false,
        }
    };
}

/// Reference symbol of each unit kind that has one. Feed, product, the
/// energy boundaries, the splitter and the logical blocks have none and
/// keep project drawings ([`fallback_code`]).
pub const TRACED: [TracedSymbol; 48] = [
    traced!(Mixer, "s4-x2673-in-line-static-mixer",
        pipes: ["path3552", "path3554"], stubs: []),
    // Heater and cooler: the process line crosses the exchanger; the
    // heating or cooling medium nozzles above and below stay drawn.
    traced!(Heater, "s2-x8079-heat-exchanger-general-condenser",
        pipes: ["path4147", "path1781"], stubs: ["path1783", "path1785"]),
    traced!(Cooler, "s2-x8079-heat-exchanger-general-condenser",
        pipes: ["path4147", "path1781"], stubs: ["path1783", "path1785"]),
    traced!(Valve, "s6-2101-valve-general",
        pipes: ["path7830", "path7832"], stubs: []),
    traced!(HeatExchanger, "s2-2511-heat-exchanger-with-straight-tubes-fixed-tube-plates",
        pipes: ["path1811", "path1813", "path1815", "path1817"], stubs: []),
    traced!(Pump, "s4-2301-pump-liquid-type-general",
        pipes: ["path3554-4", "path4086"], stubs: []),
    traced!(Compressor, "s4-2302-compressor-vacuum-pump-general",
        pipes: ["path4251", "path4253"], stubs: []),
    traced!(Expander, "s5-2571-turbine-general",
        pipes: ["path2892", "path3544"], stubs: []),
    traced!(Separator, "s1-2062-tank-vessel-with-dished-ends",
        pipes: ["path904-9", "path906-3"], stubs: []),
    // A17 palette (docs/a17-architecture.md). A symbol does not imply a
    // calculation; see `UnitKind::calculation_status`.
    traced!(ComponentSplitter, "s3-x2616-impact-separator"),
    traced!(Selector, "s6-2103-valve-three-way-type-general"),
    traced!(StreamCutter, "s6-516-reducer"),
    traced!(
        ThreePhaseSeparator,
        "s1-x8201-vessel-with-two-different-diameters"
    ),
    traced!(Decanter, "s3-x8031-gravity-separator-settling-chamber"),
    traced!(Tank, "s1-x2063-tank-with-conical-roof-and-flat-bottom"),
    traced!(
        AirCooler,
        "s2-x2505-heat-exchanger-with-finned-tube-and-fan"
    ),
    traced!(FiredHeater, "s2-2533-furnace"),
    traced!(MultiStreamExchanger, "s2-2516-heat-exchanger-of-plate-type"),
    traced!(
        MultistageCompressor,
        "s4-x8102-compressor-vacuum-pump-turbo-type"
    ),
    traced!(PipeSegment, "s6-x322-pipeline-with-thermal-insulation"),
    traced!(
        ReliefValve,
        "s6-x2124-safety-valve-spring-loaded-globe-type"
    ),
    traced!(ShortcutColumn, "s1-x8100-column-general"),
    traced!(DistillationColumn, "s1-x8101-tray-column-general"),
    traced!(
        Absorber,
        "s1-x8016-column-with-two-fixed-bed-sections-and-intermediate-spray-nozzle"
    ),
    traced!(ThreePhaseColumn, "s1-x8012-column-with-valve-trays"),
    traced!(LiquidExtractor, "s1-x8013-column-staggered-baffle-trays"),
    traced!(
        PetroleumFractionator,
        "s1-x8011-column-with-bubble-cap-trays"
    ),
    // ISO 10628-2 has no general reactor symbol. A fixed-bed column stands
    // for the stoichiometric (catalytic) reactor and a fluidized-bed column
    // for the yield reactor (typical of cracking yield shifts). The four
    // sheet pipe lines are left out; the side walls carry the nozzles.
    traced!(ConversionReactor, "s1-x8015-column-with-fixed-bed",
        pipes: ["path4628", "path4630", "path4636", "path4638"], stubs: []),
    traced!(YieldReactor, "s1-x8014-column-with-fluidized-bed",
        pipes: ["path4392", "path4394", "path4406", "path4408"], stubs: []),
    // Equilibrium reactor: vessel with heating/cooling jacket; Gibbs reactor:
    // spherical vessel; CSTR: jacketed vessel with agitator; plug-flow
    // reactor: jacketed pipeline (tube with coolant); batch reactor: vessel
    // with a heating/cooling coil.
    traced!(EquilibriumReactor, "s1-x2069-vessel-with-dished-ends-and-heating-cooling-jacket",
        pipes: ["path2229", "path2231", "path2243", "path2245"], stubs: []),
    traced!(GibbsReactor, "s1-2063-spherical-vessel",
        pipes: ["path1041", "path1049", "path1051", "path1053"], stubs: []),
    traced!(Cstr, "s1-x8006-jacketed-vessel-with-dished-ends-and-agitator-driven-by-electric-motor",
        pipes: ["path4191", "path4193", "path4199", "path4272"], stubs: []),
    traced!(PlugFlowReactor, "s6-x409-jacked-pipeline",
        pipes: ["path6507", "path6509"], stubs: []),
    traced!(BatchReactor, "s1-x8000-vessel-with-full-tube-heating-cooling-coil",
        pipes: ["path7578", "path7580", "path4115", "path4117"], stubs: []),
    traced!(Cyclone, "s3-x2618-separator-cyclone-type"),
    traced!(Hydrocyclone, "s3-x2618-separator-cyclone-type"),
    traced!(Baghouse, "s3-x8022-gas-filter-bag-candle-or-cartridge-type"),
    traced!(
        RotaryFilter,
        "s3-x8119-liquid-rotary-filter-drum-or-disc-type"
    ),
    traced!(SolidsSeparator, "s3-x8081-separator-sifter-general"),
    traced!(
        Crystallizer,
        "s1-x8008-vessel-with-dished-roof-and-conical-bottom"
    ),
    traced!(Centrifuge, "s3-x8035-centrifuge-with-solid-shell"),
    traced!(
        ElectrostaticPrecipitator,
        "s3-x8125-separator-precipitator-electrostatic-type"
    ),
    traced!(VenturiScrubber, "s3-x8034-separator-scrubber-venturi-type"),
    traced!(Crusher, "s4-x8085-crusher-general"),
    traced!(Screen, "s3-x8123-screening-device-sieve-strainer-general"),
    traced!(Dryer, "s4-c0046-drier-general"),
    traced!(Granulator, "s4-x8057-pelletizing-disc"),
    traced!(SolidsWasher, "s3-x8124-separator-scrubber-wet-type"),
];

/// Short code drawn in the box of a kind without a reference symbol:
/// ISO 10628-2 has no symbols for logical blocks. `None` for kinds with a
/// dedicated project drawing (feed, product, energy boundaries, splitter).
pub fn fallback_code(kind: UnitKind) -> Option<&'static str> {
    Some(match kind {
        UnitKind::Multiplier => "×n",
        UnitKind::Adjust => "ADJ",
        UnitKind::Set => "SET",
        UnitKind::Spreadsheet => "f(x)",
        UnitKind::Balance => "BAL",
        UnitKind::RecycleBlock => "RCY",
        UnitKind::ControlBlock => "PID",
        UnitKind::UserUnit => "USER",
        _ => return None,
    })
}

/// Stroke width of equipment outlines on the reference sheets (mm). Line
/// weights are drawn relative to it.
const BODY_WIDTH: f32 = 0.5;

/// Thinnest line drawn on screen (points), so thin reference lines stay
/// visible when zoomed out.
const MIN_WIDTH: f32 = 1.0;

/// Segments per Bézier curve when flattening.
const CURVE_STEPS: usize = 12;

/// One polyline of a traced symbol, in coordinates normalised to the
/// symbol's bounding box (0..1 on both axes).
#[derive(Debug, Clone)]
struct Line {
    points: Vec<Pos2>,
    closed: bool,
    /// Stroke width relative to [`BODY_WIDTH`].
    weight: f32,
    /// Nozzle stub outside the body.
    stub: bool,
    /// Dash and gap length as fractions of the symbol width, for lines
    /// drawn dashed on the reference sheet (fixed-bed boundaries).
    dash: Option<(f32, f32)>,
}

/// Parsed, normalised geometry of a [`TracedSymbol`].
#[derive(Debug, Clone)]
pub struct Outline {
    lines: Vec<Line>,
    /// Convex hull of the body, filled with the background colour.
    mask: Vec<Pos2>,
    /// Width over height of the bounding box.
    pub aspect: f32,
}

impl Outline {
    /// Polylines mapped into `r` (already fitted to the aspect ratio).
    fn mapped(&self, r: Rect) -> impl Iterator<Item = (&Line, Vec<Pos2>)> {
        self.lines
            .iter()
            .map(move |l| (l, l.points.iter().map(|&p| map(r, p)).collect()))
    }
}

fn map(r: Rect, p: Pos2) -> Pos2 {
    pos2(r.left() + p.x * r.width(), r.top() + p.y * r.height())
}

/// Traced symbol of `kind`, if it has one.
pub fn traced_symbol(kind: UnitKind) -> Option<&'static TracedSymbol> {
    TRACED.iter().find(|t| t.kind == kind)
}

/// Parsed outline of `kind`'s traced symbol. The SVGs are parsed once; the
/// embedded files are checked by the tests, so a parse failure here means a
/// broken build and falls back to no drawing rather than a panic.
pub fn outline(kind: UnitKind) -> Option<&'static Outline> {
    static OUTLINES: OnceLock<Vec<Option<Outline>>> = OnceLock::new();
    let all = OUTLINES.get_or_init(|| TRACED.iter().map(|t| parse(t).ok()).collect());
    let index = TRACED.iter().position(|t| t.kind == kind)?;
    all[index].as_ref()
}

/// Parses a traced symbol into normalised polylines.
pub fn parse(t: &TracedSymbol) -> Result<Outline, String> {
    let tree = usvg::Tree::from_str(t.svg, &usvg::Options::default())
        .map_err(|e| format!("{}: {e}", t.asset))?;
    let mut raw: Vec<(String, Line)> = Vec::new();
    collect(tree.root(), &mut raw);
    for id in t.pipes.iter().chain(t.stubs) {
        if !raw.iter().any(|(rid, _)| rid == id) {
            return Err(format!("{}: element {id} not found", t.asset));
        }
    }
    let mut lines: Vec<Line> = raw
        .into_iter()
        .filter(|(id, _)| !t.pipes.contains(&id.as_str()))
        .map(|(id, mut line)| {
            line.stub = t.stubs.contains(&id.as_str());
            line
        })
        .collect();
    let points = || lines.iter().flat_map(|l| l.points.iter().copied());
    let bounds = Rect::from_points(&points().collect::<Vec<_>>());
    if lines.is_empty() || bounds.width() <= 0.0 || bounds.height() <= 0.0 {
        return Err(format!("{}: no drawable geometry", t.asset));
    }
    for line in &mut lines {
        if let Some((d, g)) = &mut line.dash {
            *d /= bounds.width();
            *g /= bounds.width();
        }
        for p in &mut line.points {
            *p = pos2(
                (p.x - bounds.left()) / bounds.width(),
                (p.y - bounds.top()) / bounds.height(),
            );
        }
    }
    let body: Vec<Pos2> = lines
        .iter()
        .filter(|l| !l.stub)
        .flat_map(|l| l.points.iter().copied())
        .collect();
    Ok(Outline {
        mask: convex_hull(body),
        aspect: bounds.width() / bounds.height(),
        lines,
    })
}

/// Collects the stroked paths below `group` in drawing coordinates.
fn collect(group: &usvg::Group, out: &mut Vec<(String, Line)>) {
    for node in group.children() {
        match node {
            usvg::Node::Group(g) => collect(g, out),
            usvg::Node::Path(p) => {
                let Some(stroke) = p.stroke() else { continue };
                let Some(data) = p.data().clone().transform(p.abs_transform()) else {
                    continue;
                };
                let weight = stroke.width().get() / BODY_WIDTH;
                // Drawing units; normalised with the points in `parse`.
                let dash = stroke
                    .dasharray()
                    .filter(|d| d.len() >= 2 && d[0] > 0.0)
                    .map(|d| (d[0], d[1]));
                for (points, closed) in flatten(&data) {
                    let line = Line {
                        points,
                        closed,
                        weight,
                        stub: false,
                        dash,
                    };
                    out.push((p.id().to_owned(), line));
                }
            }
            _ => {}
        }
    }
}

/// Flattens a path into polylines, with a flag for closed subpaths.
fn flatten(data: &usvg::tiny_skia_path::Path) -> Vec<(Vec<Pos2>, bool)> {
    let pt = |p: Point| pos2(p.x, p.y);
    let mut out = Vec::new();
    let mut cur: Vec<Pos2> = Vec::new();
    let mut start = Pos2::ZERO;
    let finish = |cur: &mut Vec<Pos2>, closed: bool, out: &mut Vec<(Vec<Pos2>, bool)>| {
        // A subpath that ends where it starts is closed even without `Z`
        // (the exchanger shell, for example), so its corner is joined.
        let ends_at_start = cur.len() > 2
            && cur
                .first()
                .zip(cur.last())
                .is_some_and(|(a, b)| a.distance(*b) < 1e-3);
        if ends_at_start {
            cur.pop();
        }
        let closed = closed || ends_at_start;
        if cur.len() > 1 {
            out.push((std::mem::take(cur), closed));
        }
        cur.clear();
    };
    for seg in data.segments() {
        let last = cur.last().copied().unwrap_or(start);
        if cur.is_empty() && !matches!(seg, PathSegment::MoveTo(_)) {
            cur.push(start);
        }
        match seg {
            PathSegment::MoveTo(p) => {
                finish(&mut cur, false, &mut out);
                start = pt(p);
                cur.push(start);
            }
            PathSegment::LineTo(p) => cur.push(pt(p)),
            PathSegment::QuadTo(c, p) => {
                let (c, p) = (pt(c), pt(p));
                for i in 1..=CURVE_STEPS {
                    let t = i as f32 / CURVE_STEPS as f32;
                    let u = 1.0 - t;
                    cur.push(
                        (last.to_vec2() * u * u + c.to_vec2() * 2.0 * u * t + p.to_vec2() * t * t)
                            .to_pos2(),
                    );
                }
            }
            PathSegment::CubicTo(c1, c2, p) => {
                let (c1, c2, p) = (pt(c1), pt(c2), pt(p));
                for i in 1..=CURVE_STEPS {
                    let t = i as f32 / CURVE_STEPS as f32;
                    let u = 1.0 - t;
                    cur.push(
                        (last.to_vec2() * u * u * u
                            + c1.to_vec2() * 3.0 * u * u * t
                            + c2.to_vec2() * 3.0 * u * t * t
                            + p.to_vec2() * t * t * t)
                            .to_pos2(),
                    );
                }
            }
            PathSegment::Close => finish(&mut cur, true, &mut out),
        }
    }
    finish(&mut cur, false, &mut out);
    out
}

/// Convex hull (Andrew's monotone chain).
fn convex_hull(mut pts: Vec<Pos2>) -> Vec<Pos2> {
    pts.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    pts.dedup_by(|a, b| a.distance(*b) < 1e-6);
    if pts.len() < 3 {
        return pts;
    }
    let cross = |o: Pos2, a: Pos2, b: Pos2| (a - o).x * (b - o).y - (a - o).y * (b - o).x;
    let half = |iter: &mut dyn Iterator<Item = Pos2>| {
        let mut chain: Vec<Pos2> = Vec::new();
        for p in iter {
            while chain.len() >= 2
                && cross(chain[chain.len() - 2], chain[chain.len() - 1], p) <= 0.0
            {
                chain.pop();
            }
            chain.push(p);
        }
        chain.pop();
        chain
    };
    let mut hull = half(&mut pts.iter().copied());
    hull.extend(half(&mut pts.iter().rev().copied()));
    hull
}

/// Rectangle of aspect ratio `aspect` centred in and fitting `r`.
fn fit(r: Rect, aspect: f32) -> Rect {
    let size = if r.width() > r.height() * aspect {
        vec2(r.height() * aspect, r.height())
    } else {
        vec2(r.width(), r.width() / aspect)
    };
    Rect::from_center_size(r.center(), size)
}

/// Draws the symbol of `kind` in the screen rectangle `r`. `fill` masks
/// grid and lines behind the symbol (use transparent for palette icons).
pub fn draw_unit(painter: &Painter, kind: UnitKind, r: Rect, stroke: Stroke, fill: Color32) {
    if let Some(outline) = outline(kind) {
        draw_traced(painter, kind, outline, r, stroke, fill);
        return;
    }
    let c = r.center();
    match kind {
        UnitKind::Feed | UnitKind::Product => {
            // Project symbol: flow boundary flag pointing in the flow
            // direction. The reference has no feed/product symbol.
            let tip = r.height() * 0.5;
            let pts = vec![
                r.left_top(),
                pos2(r.right() - tip, r.top()),
                pos2(r.right(), c.y),
                pos2(r.right() - tip, r.bottom()),
                r.left_bottom(),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
        }
        UnitKind::EnergySource | UnitKind::EnergySink => {
            // Project symbol (#171): a dashed boundary flag marked Q, as the
            // dashed energy streams; the reference has no energy symbol.
            let tip = r.height() * 0.5;
            let pts = [
                r.left_top(),
                pos2(r.right() - tip, r.top()),
                pos2(r.right(), c.y),
                pos2(r.right() - tip, r.bottom()),
                r.left_bottom(),
            ];
            painter.add(Shape::convex_polygon(pts.to_vec(), fill, Stroke::NONE));
            let mut closed = pts.to_vec();
            closed.push(pts[0]);
            let dash = (r.height() * 0.2).max(2.0);
            painter.extend(Shape::dashed_line(&closed, stroke, dash, dash * 0.6));
            painter.text(
                pos2(c.x - tip * 0.3, c.y),
                eframe::egui::Align2::CENTER_CENTER,
                "Q",
                eframe::egui::FontId::monospace((r.height() * 0.55).max(6.0)),
                stroke.color,
            );
        }
        UnitKind::Splitter => {
            // Project symbol: one inlet forking to the two outlet nozzles.
            // ISO 10628-2 has no unit symbol for a stream split.
            painter.rect(r, 0, fill, stroke, eframe::egui::StrokeKind::Middle);
            let thin = Stroke::new(stroke.width * 0.7, stroke.color);
            let y1 = r.top() + r.height() / 3.0;
            let y2 = r.top() + r.height() * 2.0 / 3.0;
            let fork = pos2(c.x, c.y);
            painter.line_segment([pos2(r.left(), c.y), fork], thin);
            painter.line(
                vec![pos2(r.right(), y1), pos2(c.x + r.width() * 0.2, y1), fork],
                thin,
            );
            painter.line(
                vec![pos2(r.right(), y2), pos2(c.x + r.width() * 0.2, y2), fork],
                thin,
            );
        }
        // Logical blocks: a box with a short code. Every other kind has a
        // traced symbol; an outline that failed to parse (caught by the
        // tests) falls back to a plain box.
        _ => {
            painter.rect(r, 0, fill, stroke, eframe::egui::StrokeKind::Middle);
            if let Some(code) = fallback_code(kind) {
                let size = (r.height() * 0.38).min(r.width() / code.len().max(2) as f32 * 1.4);
                painter.text(
                    c,
                    eframe::egui::Align2::CENTER_CENTER,
                    code,
                    eframe::egui::FontId::monospace(size.max(5.0)),
                    stroke.color,
                );
            }
        }
    }
}

fn draw_traced(
    painter: &Painter,
    kind: UnitKind,
    outline: &Outline,
    r: Rect,
    stroke: Stroke,
    fill: Color32,
) {
    let r = fit(r, outline.aspect);
    if fill.a() > 0 && outline.mask.len() >= 3 {
        let mask = outline.mask.iter().map(|&p| map(r, p)).collect();
        painter.add(Shape::convex_polygon(mask, fill, Stroke::NONE));
    }
    let mut stubs: Vec<Vec<Pos2>> = Vec::new();
    for (line, points) in outline.mapped(r) {
        let width = (stroke.width * line.weight).max(stroke.width.min(MIN_WIDTH));
        let s = Stroke::new(width, stroke.color);
        if line.stub {
            stubs.push(points.clone());
        }
        // A zero-length stroke with a round cap is a dot in SVG (the
        // fluidized-bed particles); egui draws nothing for it.
        if points.windows(2).all(|w| w[0].distance(w[1]) < 1e-3) {
            painter.circle_filled(points[0], width / 2.0, stroke.color);
            continue;
        }
        match line.dash {
            Some((d, g)) if !line.closed => {
                painter.extend(Shape::dashed_line(&points, s, d * r.width(), g * r.width()))
            }
            _ => {
                painter.add(if line.closed {
                    Shape::closed_line(points, s)
                } else {
                    Shape::line(points, s)
                });
            }
        }
    }
    // Project addition: an arrow on the medium shows heat in (heater, into
    // the exchanger through the top nozzle) or out (cooler, away from it
    // through the bottom nozzle).
    if matches!(kind, UnitKind::Heater | UnitKind::Cooler) {
        let lowest = |l: &Vec<Pos2>| l.iter().map(|p| p.y).fold(f32::MIN, f32::max);
        let stub = if kind == UnitKind::Heater {
            stubs.iter().min_by(|a, b| lowest(a).total_cmp(&lowest(b)))
        } else {
            stubs.iter().max_by(|a, b| lowest(a).total_cmp(&lowest(b)))
        };
        if let Some(tip) = stub.and_then(|l| l.iter().copied().max_by(|a, b| a.y.total_cmp(&b.y))) {
            let head = (r.width() / 2.0 * 0.35).max(4.0);
            arrowhead(painter, tip, vec2(0.0, 1.0), head, stroke.color);
        }
    }
}

/// Scales a symbol to fit an icon area, keeping its aspect ratio.
pub fn icon_rect(area: Rect, kind: UnitKind) -> Rect {
    let size = crate::view::block_size(kind);
    let scale = (area.width() / size.x).min(area.height() / size.y);
    Rect::from_center_size(area.center(), size * scale)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::view::{block_rect, block_size, distance_to_polyline, port_pos};

    const MANIFEST: &str = include_str!("../../../docs/design-reference/iso-10628/manifest.json");

    #[test]
    fn every_traced_symbol_parses_and_is_in_the_manifest() {
        for t in &TRACED {
            let o = parse(t).unwrap_or_else(|e| panic!("{e}"));
            assert!(outline(t.kind).is_some(), "{}", t.kind);
            assert!(
                MANIFEST.contains(&format!("\"id\": \"{}\"", t.asset)),
                "{} not in manifest",
                t.asset
            );
            assert!(
                t.svg.contains("CC BY-SA 4.0"),
                "{} lost its license",
                t.asset
            );
            assert!(!o.lines.is_empty(), "{}: no lines", t.asset);
            assert!(o.mask.len() >= 3, "{}: no mask", t.asset);
        }
    }

    #[test]
    fn unit_kinds_without_a_reference_symbol_are_the_project_ones() {
        let missing: Vec<UnitKind> = UnitKind::ALL
            .into_iter()
            .filter(|&k| traced_symbol(k).is_none())
            .collect();
        // Logical blocks have no ISO 10628 symbol; everything else does.
        assert_eq!(
            missing,
            [
                UnitKind::Feed,
                UnitKind::Product,
                UnitKind::EnergySource,
                UnitKind::EnergySink,
                UnitKind::Splitter,
                UnitKind::Multiplier,
                UnitKind::Adjust,
                UnitKind::Set,
                UnitKind::Spreadsheet,
                UnitKind::Balance,
                UnitKind::RecycleBlock,
                UnitKind::ControlBlock,
                UnitKind::UserUnit,
            ]
        );
        for kind in missing {
            let project_drawing = matches!(
                kind,
                UnitKind::Feed
                    | UnitKind::Product
                    | UnitKind::EnergySource
                    | UnitKind::EnergySink
                    | UnitKind::Splitter
            );
            assert_eq!(fallback_code(kind).is_none(), project_drawing, "{kind}");
        }
    }

    /// The block keeps the reference proportions, so the traced geometry
    /// fills it exactly when scaled uniformly.
    #[test]
    fn blocks_have_the_proportions_of_their_symbol() {
        for t in &TRACED {
            let o = outline(t.kind).unwrap();
            let size = block_size(t.kind);
            let block = size.x / size.y;
            assert!(
                (block / o.aspect - 1.0).abs() < 0.01,
                "{}: block {block}, symbol {}",
                t.kind,
                o.aspect
            );
        }
    }

    /// Ports are defined in `view.rs`, apart from the drawings; each one
    /// must still land on the traced line work.
    #[test]
    fn every_port_lies_on_the_traced_line_work() {
        for t in TRACED.iter().filter(|t| t.reviewed_ports) {
            let o = outline(t.kind).unwrap();
            let rect = block_rect(t.kind, Pos2::ZERO);
            let lines: Vec<Vec<Pos2>> = o
                .mapped(fit(rect, o.aspect))
                .filter(|(l, _)| !l.stub)
                .map(|(l, mut pts)| {
                    if l.closed {
                        pts.push(pts[0]);
                    }
                    pts
                })
                .collect();
            for i in 0..t.kind.ports().len() {
                let p = port_pos(t.kind, Pos2::ZERO, i);
                let d = lines
                    .iter()
                    .map(|l| distance_to_polyline(p, l))
                    .fold(f32::INFINITY, f32::min);
                assert!(
                    d < 0.5,
                    "{} port {i} at {p:?} is {d} from the symbol",
                    t.kind
                );
            }
        }
    }

    #[test]
    fn mask_covers_the_body_and_stays_inside_the_block() {
        for t in &TRACED {
            let o = outline(t.kind).unwrap();
            for p in &o.mask {
                assert!(
                    (-1e-4..=1.0 + 1e-4).contains(&p.x) && (-1e-4..=1.0 + 1e-4).contains(&p.y),
                    "{}: {p:?}",
                    t.asset
                );
            }
        }
    }

    #[test]
    fn pipe_lines_from_the_sheet_are_left_out() {
        let t = traced_symbol(UnitKind::Pump).unwrap();
        let o = parse(t).unwrap();
        // Circle and discharge chevron only; the reference's two pipe
        // stubs are replaced by the streams.
        assert_eq!(o.lines.len(), 2);
        assert!(o.lines.iter().all(|l| (l.weight - 1.0).abs() < 1e-6));
    }
}
