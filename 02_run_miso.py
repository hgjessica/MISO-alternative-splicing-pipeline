#!/usr/bin/env python3
"""Step 2: compute Psi per sample and event type with `miso --run`.

    python 02_run_miso.py --config config.json run --sample sample_1 [--event-type SE]
    python 02_run_miso.py --config config.json insert-len --sample sample_1
    python 02_run_miso.py --config config.json submit run

run: miso --run <index/TYPE> <bam> --output-dir miso/<TYPE>/<sample> --read-len N, plus the optional
flags from config "miso" (paired_end, overhang_len, num_processes, settings_filename,
no_filter_events, prefilter, pass_event_type). The read length comes from sample.read_len, else
miso.read_len. Under SLURM one array task = one (sample, event type) pair.
insert-len (optional, paired-end only): estimate the insert-size mean/sd with MISO's
exon_utils + pe_utils; copy the values into miso.paired_end. Output: insert_len/<sample>/
Output of run: miso/<TYPE>/<sample>/<chrom>/<event>.miso
"""
import sys
from pathlib import Path

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array

STEP_UNITS = {"run": "sample_event"}


def step_run(ctx, args):
    """miso --run for one (sample, event type)."""
    raise_open_files()
    m = ctx.miso
    for s, t in pick_tasks(ctx, "sample_event", args.sample, args.event_type):
        task_tmp(ctx, "miso")
        idx, bam = ctx.index_dir(t), Path(s["bam"])
        need(ctx, idx, bam, str(bam) + ".bai")
        out = ctx.miso_dir(t, s["name"])
        if not ctx.dry:
            out.mkdir(parents=True, exist_ok=True)
        cmd = ["miso", "--run", idx, bam, "--output-dir", out, "--read-len", read_len(ctx, s)]
        pe = s.get("paired_end", m.get("paired_end"))
        if pe:
            cmd += ["--paired-end", pe["mean"], pe["sd"]]
        if m.get("overhang_len") is not None:
            cmd += ["--overhang-len", m["overhang_len"]]
        if m.get("settings_filename"):
            cmd += ["--settings-filename", m["settings_filename"]]
        if m.get("num_processes"):
            cmd += ["-p", m["num_processes"]]
        if m.get("no_filter_events"):
            cmd.append("--no-filter-events")
        if m.get("prefilter"):
            cmd.append("--prefilter")
        if m.get("pass_event_type"):
            cmd += ["--event-type", t]
        if run(ctx, cmd):
            sys.exit(f"[ERROR] miso failed for {s['name']} {t}")
        print(f"[OK] {s['name']} {t}: {out}")


def step_insert_len(ctx, args):
    """Estimate the paired-end insert length (optional, not part of the single-end default)."""
    samples = kept_samples(ctx.cfg) if args.sample == "ALL" else [sample_by_name(ctx.cfg, args.sample)]
    etype = args.event_type or ctx.types[0]
    gff = ctx.cfg["annotation"]["events"][etype]
    need(ctx, gff)
    const = ctx.path("insert_len", "const_exons")
    if not ctx.dry:
        const.mkdir(parents=True, exist_ok=True)
    if run(ctx, ["python", "-m", "misopy.exon_utils", "--get-const-exons", gff,
                 "--min-exon-size", args.min_exon_size, "--output-dir", const]):
        sys.exit("[ERROR] exon_utils failed")
    for s in samples:
        out = ctx.path("insert_len", s["name"])
        need(ctx, s["bam"])
        if not ctx.dry:
            out.mkdir(parents=True, exist_ok=True)
        cexons = const / (Path(gff).stem + ".min_" + args.min_exon_size + ".const_exons.gff")
        if run(ctx, ["python", "-m", "misopy.pe_utils", "--compute-insert-len",
                     s["bam"], cexons, "--output-dir", out]):
            sys.exit(f"[ERROR] pe_utils failed for {s['name']}")
        print(f"[OK] {s['name']}: see the '#mean=..,sdev=..' line in {out}/*.insert_len")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="miso --run")
    add_selectors(p, "sample")
    p.set_defaults(fn=step_run)
    p = sub.add_parser("insert-len", help="estimate insert length (paired-end)")
    p.add_argument("--sample", required=True, help="sample name or ALL")
    p.add_argument("--event-type", help="event GFF3 to take constitutive exons from (default: first)")
    p.add_argument("--min-exon-size", default="500")
    p.set_defaults(fn=step_insert_len)
    add_submit(sub, STEP_UNITS).set_defaults(
        fn=lambda c, a: submit_array(c, a, __file__, STEP_UNITS))
    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
