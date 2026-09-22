//! Markdown summary of one or more results logs: what's excluded, what was found,
//! and what's still unresolved.

use crate::log::{Kind, Record};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

/// (base, period, degrees, extras) — one formula shape, independent of C and precision.
type ShapeKey = (u32, u32, [u32; 2], String);

fn key(r: &Record) -> ShapeKey {
    let p = &r.params;
    (p.base, p.period, p.degrees, extras(r))
}

fn extras(r: &Record) -> String {
    if r.params.extras.is_empty() {
        "-".into()
    } else {
        r.params.extras.join(",")
    }
}

/// One-line shape label, e.g. `b=10 m=12 s=1..2 x=log2,log5`.
pub fn describe(r: &Record) -> String {
    let p = &r.params;
    format!(
        "b={} m={} s={}..{} x={}",
        p.base,
        p.period,
        p.degrees[0],
        p.degrees[1],
        extras(r)
    )
}

/// `c1·col1 + c2·col2 + ... = 0` over the nonzero coefficients.
pub fn formula(columns: &[String], relation: &[String]) -> String {
    columns
        .iter()
        .zip(relation)
        .filter(|(_, c)| c.as_str() != "0")
        .map(|(n, c)| format!("{c}·{n}"))
        .collect::<Vec<_>>()
        .join(" + ")
        + " = 0"
}

fn kind_name(k: Kind) -> String {
    serde_json::to_value(k)
        .ok()
        .and_then(|v| v.as_str().map(String::from))
        .unwrap_or_default()
}

pub fn render(records: &[Record]) -> String {
    let mut out = String::from("# pihunt exclusion report\n\n");
    if records.is_empty() {
        out.push_str("No records.\n");
        return out;
    }

    let mut by_shape: BTreeMap<ShapeKey, Vec<&Record>> = BTreeMap::new();
    for r in records {
        by_shape.entry(key(r)).or_default().push(r);
    }

    // Strongest exclusion per shape: largest C, then fewest dropped columns.
    let exclusions: Vec<&Record> = by_shape
        .values()
        .filter_map(|rs| {
            rs.iter()
                .filter(|r| r.outcome == Kind::Excluded)
                .max_by_key(|r| (r.params.coeff_bound, std::cmp::Reverse(r.dropped.len())))
                .copied()
        })
        .collect();

    let mut hits: BTreeMap<(ShapeKey, String), &Record> = BTreeMap::new();
    for r in records.iter().filter(|r| r.outcome == Kind::Hit) {
        let f = formula(&r.columns, r.relation.as_deref().unwrap_or_default());
        hits.entry((key(r), f)).or_insert(r);
    }
    let new_hits = hits
        .values()
        .filter(|r| r.tag.as_deref() == Some("NEW"))
        .count();

    let mut basis: BTreeSet<(u32, u32, String, String)> = BTreeSet::new();
    for r in records {
        for d in &r.dropped {
            let tag = d.tag.clone().unwrap_or_else(|| "-".into());
            basis.insert((
                r.params.base,
                r.params.period,
                formula(&r.columns, &d.relation),
                tag,
            ));
        }
    }

    // Shapes nothing settled: report the highest-precision attempt.
    let unresolved: Vec<&Record> = by_shape
        .values()
        .filter(|rs| {
            !rs.iter()
                .any(|r| matches!(r.outcome, Kind::Excluded | Kind::Hit))
        })
        .filter_map(|rs| rs.iter().max_by_key(|r| r.params.precision_digits).copied())
        .collect();

    let finders: BTreeSet<&str> = records.iter().map(|r| r.params.finder.as_str()).collect();
    let _ = writeln!(
        out,
        "Records: {} · Shapes searched: {} · Shapes excluded: {} · Hits: {} ({} NEW) · Unresolved: {} · Finders: {}\n",
        records.len(),
        by_shape.len(),
        exclusions.len(),
        hits.len(),
        new_hits,
        unresolved.len(),
        finders.into_iter().collect::<Vec<_>>().join(", ")
    );
    out.push_str(
        "An exclusion means PSLQ proved no integer relation between π and the shape's reduced basis \
         has every |coefficient| ≤ C. Columns listed as dropped were removed first as rational \
         combinations of the others (see Basis relations).\n\n",
    );

    if new_hits > 0 {
        out.push_str(
            "## 🚨 NEW HITS\n\nVerified at 2× precision and not in the known-formula table.\n\n",
        );
        for r in hits.values().filter(|r| r.tag.as_deref() == Some("NEW")) {
            let f = formula(&r.columns, r.relation.as_deref().unwrap_or_default());
            let _ = writeln!(out, "- **{}**: `{f}` (job `{}`)", describe(r), r.job_id);
        }
        out.push('\n');
    }

    if !hits.is_empty() {
        out.push_str("## Hits\n\n| base | period | degrees | extras | tag | relation |\n|---|---|---|---|---|---|\n");
        for ((k, f), r) in &hits {
            let _ = writeln!(
                out,
                "| {} | {} | {}..{} | {} | {} | `{f}` |",
                k.0,
                k.1,
                k.2[0],
                k.2[1],
                k.3,
                r.tag.as_deref().unwrap_or("-")
            );
        }
        out.push('\n');
    }

    out.push_str("## Exclusions\n\n");
    if exclusions.is_empty() {
        out.push_str("None yet.\n\n");
    } else {
        out.push_str(
            "| base | period | degrees | extras | max \\|coeff\\| ≤ | digits | dropped | finder |\n\
             |---|---|---|---|---|---|---|---|\n",
        );
        for r in &exclusions {
            let (b, m, d, x) = key(r);
            let dropped = if r.dropped.is_empty() {
                "-".to_string()
            } else {
                r.dropped
                    .iter()
                    .map(|d| d.column.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            let _ = writeln!(
                out,
                "| {b} | {m} | {}..{} | {x} | {} | {} | {dropped} | {} |",
                d[0], d[1], r.params.coeff_bound, r.params.precision_digits, r.params.finder
            );
        }
        out.push('\n');
    }

    if !basis.is_empty() {
        out.push_str(
            "## Basis relations\n\n| base | period | relation | tag |\n|---|---|---|---|\n",
        );
        for (b, m, f, tag) in &basis {
            let _ = writeln!(out, "| {b} | {m} | `{f}` | {tag} |");
        }
        out.push('\n');
    }

    out.push_str("## Unresolved\n\n");
    if unresolved.is_empty() {
        out.push_str("None.\n");
    } else {
        out.push_str(
            "| base | period | degrees | extras | outcome | digits | note |\n|---|---|---|---|---|---|---|\n",
        );
        for r in &unresolved {
            let (b, m, d, x) = key(r);
            let _ = writeln!(
                out,
                "| {b} | {m} | {}..{} | {x} | {} | {} | {} |",
                d[0],
                d[1],
                kind_name(r.outcome),
                r.params.precision_digits,
                r.note.as_deref().unwrap_or("-")
            );
        }
    }
    out
}
