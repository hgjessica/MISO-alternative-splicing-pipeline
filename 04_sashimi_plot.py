#!/usr/bin/env python3
"""Step 4 (optional): sashimi plots for selected events and a Bayes factor distribution plot.

    python 04_sashimi_plot.py --config config.json plot --comparison treated_vs_control [--event-type SE]
    python 04_sashimi_plot.py --config config.json bf-dist --comparison treated_vs_control
    python 04_sashimi_plot.py --config config.json submit plot

plot: writes a MISO plot-settings file (control first, treated second) from config "plots", then
runs sashimi_plot --plot-event <event> index/<TYPE> <settings> for each event: those listed in
plots.events (matching the event type) plus the top plots.top_n significant events by Bayes
factor from step 3. bf-dist: sashimi_plot --plot-bf-dist <.miso_bf> <settings>.
Output: sashimi/<TYPE>/<comparison>/<event>.pdf, sashimi/<comparison>.settings.txt
"""
import csv
import subprocess
import sys
from pathlib import Path

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array

STEP_UNITS = {"plot": "comparison_event"}

# plot parameters taken from the original MISO settings files; override in config "plots"
PLOT_DEFAULTS = {
    "fig_width": 7, "fig_height": 5, "intron_scale": 30, "exon_scale": 4, "logged": False,
    "font_size": 6, "ymax": 150, "show_posteriors": True, "bar_posteriors": False,
    "number_junctions": True, "resolution": ".5", "posterior_bins": 40,
    "gene_posterior_ratio": 5, "colors": ["#CC0011", "#FF8800"], "bar_color": "b",
    "bf_thresholds": [0, 1, 2, 5, 10, 20],
}


def mapped_reads(ctx, sample):
    """Mapped-read count used by MISO to normalise read density (primary mapped alignments;
    override with sample.mapped_reads)."""
    if sample.get("mapped_reads"):
        return int(sample["mapped_reads"])
    if ctx.dry:
        return 0
    out = subprocess.run(["samtools", "view", "-c", "-F", "260", sample["bam"]],
                         capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"[ERROR] samtools view failed for {sample['bam']}")
    return int(out.stdout.strip())


def write_settings(ctx, comp, etype):
    """Write the MISO plotting settings file for one comparison and event type."""
    p = dict(PLOT_DEFAULTS, **ctx.cfg.get("plots", {}).get("settings", {}))
    ctrl, trt = (sample_by_name(ctx.cfg, comp[k]) for k in ("control", "treated"))
    q = lambda xs: "[" + ", ".join(f'"{x}"' for x in xs) + "]"
    lines = ["[data]", 'bam_prefix = /', f"bam_files = {q([ctrl['bam'], trt['bam']])}",
             f"miso_prefix = {ctx.path('miso', etype)}/",
             f"miso_files = {q([ctrl['name'], trt['name']])}", "", "[plotting]"]
    for k in ("fig_width", "fig_height", "intron_scale", "exon_scale", "logged", "font_size",
              "ymax", "show_posteriors", "bar_posteriors", "number_junctions", "resolution",
              "posterior_bins", "gene_posterior_ratio"):
        lines.append(f"{k} = {p[k]}")
    lines += [f"colors = {q(p['colors'])}",
              f"coverages = [{mapped_reads(ctx, ctrl)}, {mapped_reads(ctx, trt)}]",
              f'bar_color = "{p["bar_color"]}"', f"bf_thresholds = {p['bf_thresholds']}", ""]
    out = ctx.path("sashimi", f"{comp['name']}.{etype}.settings.txt")
    if not ctx.dry:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(lines))
    return out


def events_to_plot(ctx, comp, etype):
    plots = ctx.cfg.get("plots", {})
    events = [e["event"] for e in plots.get("events", []) if e.get("type") == etype]
    sig = ctx.path("filtered", etype, f"{comp['name']}.significant.tsv")
    n = int(plots.get("top_n", 0))
    if n and sig.exists():
        with open(sig, newline="") as fh:
            events += [r["event_name"] for _, r in zip(range(n), csv.DictReader(fh, delimiter="\t"))]
    elif n and not ctx.dry:
        print(f"[WARN] {sig} not found; run step 3 filter for top_n events")
    return list(dict.fromkeys(events))


def step_plot(ctx, args):
    """sashimi_plot --plot-event for one (comparison, event type)."""
    env = {"MPLBACKEND": "Agg"}
    for c, t in pick_tasks(ctx, "comparison_event", args.comparison, args.event_type):
        task_tmp(ctx, "sashimi")
        idx = ctx.index_dir(t)
        need(ctx, idx)
        settings = write_settings(ctx, c, t)
        out = ctx.path("sashimi", t, c["name"])
        if not ctx.dry:
            out.mkdir(parents=True, exist_ok=True)
        events = events_to_plot(ctx, c, t)
        if not events:
            print(f"[SKIP] {c['name']} {t}: no events (set plots.events or plots.top_n)")
        for ev in events:
            cmd = ["sashimi_plot", "--plot-event", ev, idx, settings, "--output-dir", out]
            if run(ctx, cmd, env=env):
                print(f"[WARN] sashimi_plot failed for {ev}")


def step_bf_dist(ctx, args):
    """Bayes factor distribution plot."""
    env = {"MPLBACKEND": "Agg"}
    for c, t in pick_tasks(ctx, "comparison_event", args.comparison, args.event_type):
        bf = ctx.bf_file(t, c)
        need(ctx, bf)
        settings = write_settings(ctx, c, t)
        out = ctx.path("sashimi", t, c["name"])
        if not ctx.dry:
            out.mkdir(parents=True, exist_ok=True)
        if run(ctx, ["sashimi_plot", "--plot-bf-dist", bf, settings, "--output-dir", out], env=env):
            sys.exit(f"[ERROR] bf-dist plot failed for {c['name']} {t}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("plot", help="sashimi plots for selected events")
    add_selectors(p, "comparison")
    p.set_defaults(fn=step_plot)
    p = sub.add_parser("bf-dist", help="Bayes factor distribution plot")
    add_selectors(p, "comparison")
    p.set_defaults(fn=step_bf_dist)
    add_submit(sub, STEP_UNITS).set_defaults(
        fn=lambda c, a: submit_array(c, a, __file__, STEP_UNITS))
    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
