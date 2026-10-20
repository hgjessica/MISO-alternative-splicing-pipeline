#!/usr/bin/env python3
"""Step 3: summarize each sample, compare conditions, and keep the significant events.

    python 03_summarize_compare_filter.py --config config.json summarize --sample sample_1
    python 03_summarize_compare_filter.py --config config.json compare --comparison treated_vs_control
    python 03_summarize_compare_filter.py --config config.json filter --comparison ALL
    python 03_summarize_compare_filter.py --config config.json submit compare

summarize: summarize_miso --summarize-samples miso/<TYPE>/<sample> summary/<TYPE>/<sample>
compare:   compare_miso --compare-samples <control dir> <treated dir> compare/<TYPE>/<comparison>
           (sample 1 = control, sample 2 = treated, so diff = PSI_control - PSI_treated)
filter:    pure Python re-implementation of `filter_events` for two-isoform events, reading the
           .miso_bf file and applying the thresholds in config "filter". Also labels the
           direction of change and splits the Psi/read-count columns.
Output: summary/<TYPE>/<sample>/summary/<sample>.miso_summary
        compare/<TYPE>/<comparison>/bayes-factors/<control>_vs_<treated>.miso_bf
        filtered/<TYPE>/<comparison>.all_events.tsv (every event, with pass flag)
        filtered/<TYPE>/<comparison>.significant.tsv, filtered/filter_summary.tsv
"""
import csv
import re
import sys
from pathlib import Path

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array

STEP_UNITS = {"summarize": "sample_event", "compare": "comparison_event"}

# MISO's own recommended example thresholds (MISO docs, filter_events); the original analysis
# scripts contain no explicit threshold, so adjust to your study.
FILTER_DEFAULTS = {"delta_psi": 0.20, "bayes_factor": 10.0, "num_inc": 1, "num_exc": 1,
                   "num_sum_inc_exc": 10, "num_total": 0, "apply_both": False}

OUT_COLS = ["event_type", "comparison", "event_name", "chrom", "strand", "psi_control",
            "psi_treated", "diff", "bayes_factor", "direction",
            "control_inc", "control_exc", "control_both",
            "treated_inc", "treated_exc", "treated_both", "pass"]


def step_summarize(ctx, args):
    """summarize_miso for one (sample, event type)."""
    for s, t in pick_tasks(ctx, "sample_event", args.sample, args.event_type):
        src = ctx.miso_dir(t, s["name"])
        need(ctx, src)
        out = ctx.path("summary", t, s["name"])
        if not ctx.dry:
            out.mkdir(parents=True, exist_ok=True)
        if run(ctx, ["summarize_miso", "--summarize-samples", src, out]):
            sys.exit(f"[ERROR] summarize_miso failed for {s['name']} {t}")
        print(f"[OK] {out}/summary/{s['name']}.miso_summary")


def step_compare(ctx, args):
    """compare_miso (Bayes factors) for one (comparison, event type)."""
    raise_open_files()
    for c, t in pick_tasks(ctx, "comparison_event", args.comparison, args.event_type):
        a, b = ctx.miso_dir(t, c["control"]), ctx.miso_dir(t, c["treated"])
        need(ctx, a, b)
        out = ctx.compare_dir(t, c["name"])
        if not ctx.dry:
            out.mkdir(parents=True, exist_ok=True)
        if run(ctx, ["compare_miso", "--compare-samples", a, b, out]):
            sys.exit(f"[ERROR] compare_miso failed for {c['name']} {t}")
        print(f"[OK] {ctx.bf_file(t, c)}")


def parse_counts(text):
    """MISO count string '(0,0):99,(0,1):13,(1,0):23,(1,1):3' -> (inc, exc, both), or None
    for events with more than two isoforms. (1,0)=inclusion-only, (0,1)=exclusion-only."""
    pairs = re.findall(r"\(([01]),([01])\):(\d+)", text or "")
    if not pairs or re.search(r"\([01],[01],", text):
        return None
    d = {(a, b): int(n) for a, b, n in pairs}
    return d.get(("1", "0"), 0), d.get(("0", "1"), 0), d.get(("1", "1"), 0)


def reads_ok(cnt, f):
    inc, exc, both = cnt
    return (inc + exc + both >= f["num_total"] and inc + exc >= f["num_sum_inc_exc"]
            and inc >= f["num_inc"] and exc >= f["num_exc"])


def parse_event(row):
    """One .miso_bf row -> dict, or None if not a two-isoform event."""
    c1, c2 = parse_counts(row["sample1_counts"]), parse_counts(row["sample2_counts"])
    if c1 is None or c2 is None:
        return None
    try:
        diff, bf = float(row["diff"]), float(row["bayes_factor"].split(",")[0])
    except ValueError:
        return None
    return {"event_name": row["event_name"], "chrom": row.get("chrom", ""),
            "strand": row.get("strand", ""), "psi_control": row["sample1_posterior_mean"],
            "psi_treated": row["sample2_posterior_mean"], "diff": diff, "bayes_factor": bf,
            "direction": "control_higher" if diff > 0 else "treated_higher" if diff < 0 else "none",
            "c1": c1, "c2": c2}


def passes(ev, f):
    if abs(ev["diff"]) < f["delta_psi"] or abs(ev["bayes_factor"]) < f["bayes_factor"]:
        return False
    r1, r2 = reads_ok(ev["c1"], f), reads_ok(ev["c2"], f)
    return (r1 and r2) if f["apply_both"] else (r1 or r2)


def step_filter(ctx, args):
    """Apply the thresholds to each comparison's .miso_bf file."""
    f = dict(FILTER_DEFAULTS, **ctx.cfg.get("filter", {}))
    summary = []
    for c, t in pick_tasks(ctx, "comparison_event", args.comparison or "ALL", args.event_type):
        bf = ctx.bf_file(t, c)
        if not ctx.dry and not bf.exists():
            print(f"[WARN] missing {bf}; run compare first (skipped)")
            continue
        if ctx.dry:
            print(f"[dry-run] would filter {bf} with {f}")
            continue
        n_all = n_skip = 0
        rows_all, rows_sig = [], []
        with open(bf, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                n_all += 1
                ev = parse_event(row)
                if ev is None:
                    n_skip += 1
                    continue
                ok = passes(ev, f)
                out = dict(ev, event_type=t, comparison=c["name"], **{"pass": int(ok)})
                for pre, key in (("control", "c1"), ("treated", "c2")):
                    out[pre + "_inc"], out[pre + "_exc"], out[pre + "_both"] = ev[key]
                rows_all.append(out)
                if ok:
                    rows_sig.append(out)
        rows_sig.sort(key=lambda r: -r["bayes_factor"])
        d = ctx.path("filtered", t)
        d.mkdir(parents=True, exist_ok=True)
        for name, rows in (("all_events", rows_all), ("significant", rows_sig)):
            with open(d / f"{c['name']}.{name}.tsv", "w", newline="") as oh:
                w = csv.DictWriter(oh, OUT_COLS, delimiter="\t", extrasaction="ignore")
                w.writeheader()
                w.writerows(rows)
        up = sum(r["direction"] == "control_higher" for r in rows_sig)
        summary.append([t, c["name"], n_all, n_skip, len(rows_sig), up, len(rows_sig) - up])
        print(f"[OK] {t} {c['name']}: {len(rows_sig)}/{n_all} events pass (skipped {n_skip})")
    if summary:
        p = ctx.path("filtered", "filter_summary.tsv")
        new = not p.exists()
        with open(p, "a", newline="") as oh:
            w = csv.writer(oh, delimiter="\t")
            if new:
                w.writerow(["event_type", "comparison", "events", "skipped_multi_isoform_or_bad",
                            "significant", "control_higher", "treated_higher"])
            w.writerows(summary)


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("summarize", help="summarize_miso per sample")
    add_selectors(p, "sample")
    p.set_defaults(fn=step_summarize)
    p = sub.add_parser("compare", help="compare_miso per comparison")
    add_selectors(p, "comparison")
    p.set_defaults(fn=step_compare)
    p = sub.add_parser("filter", help="keep events passing the delta-PSI / Bayes factor filter")
    add_selectors(p, "comparison")
    p.set_defaults(fn=step_filter)
    add_submit(sub, STEP_UNITS).set_defaults(
        fn=lambda c, a: submit_array(c, a, __file__, STEP_UNITS))
    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
