#!/usr/bin/env python3
"""Step 1: check inputs and index the alternative-splicing event annotations for MISO.

    python 01_prepare_events.py --config config.json check
    python 01_prepare_events.py --config config.json index [--event-type SE]

check: every sample BAM exists, has an index, and its chromosome names match the event GFF3s.
index: index_gff --index <events GFF3> <work_dir>/index/<TYPE> for each configured event type
(SE, A3SS, A5SS, MXE, RI, ...). Run once; the index is shared by all samples.
Output: index/<TYPE>/<chrom>/*.pickle, genes_to_filenames.shelve*, genes.gff
"""
import subprocess
import sys
from pathlib import Path

from common import *  # noqa: F401,F403  (helpers shared by all steps)


def gff_chroms(path):
    seen = set()
    with open(path) as fh:
        for line in fh:
            if line and not line.startswith("#"):
                seen.add(line.split("\t", 1)[0])
    return seen


def bam_chroms(bam):
    out = subprocess.run(["samtools", "idxstats", str(bam)], capture_output=True, text=True)
    if out.returncode != 0:
        return None
    return {l.split("\t")[0] for l in out.stdout.splitlines() if l.strip()} - {"*"}


def step_check(ctx, args):
    """Report missing BAMs/indexes and chromosome-name mismatches (e.g. chr1 vs 1)."""
    bad = 0
    gff_sets = {}
    for t, gff in ctx.cfg["annotation"]["events"].items():
        if not Path(gff).exists():
            print(f"[MISSING] {t} annotation: {gff}")
            bad += 1
        else:
            gff_sets[t] = gff_chroms(gff)
    for s in kept_samples(ctx.cfg):
        bam = Path(s["bam"])
        if not bam.exists():
            print(f"[MISSING] {s['name']}: {bam}")
            bad += 1
            continue
        if not (Path(str(bam) + ".bai").exists() or bam.with_suffix(".bai").exists()):
            print(f"[NO INDEX] {s['name']}: run samtools index {bam}")
            bad += 1
            continue
        chroms = bam_chroms(bam)
        if chroms is None:
            print(f"[WARN] {s['name']}: samtools idxstats failed")
            bad += 1
            continue
        for t, g in gff_sets.items():
            if not (chroms & g):
                print(f"[MISMATCH] {s['name']}: no chromosome shared with {t} GFF3 "
                      f"(BAM e.g. {sorted(chroms)[:2]}, GFF3 e.g. {sorted(g)[:2]})")
                bad += 1
        print(f"[OK] {s['name']}: read_len {s.get('read_len', ctx.miso.get('read_len', '?'))}")
    sys.exit(1 if bad else 0)


def step_index(ctx, args):
    """index_gff for each event type."""
    types = [args.event_type] if args.event_type else ctx.types
    for t in types:
        if t not in ctx.types:
            sys.exit(f"unknown event type {t}")
        gff = ctx.cfg["annotation"]["events"][t]
        need(ctx, gff)
        out = ctx.index_dir(t)
        if not ctx.dry:
            out.mkdir(parents=True, exist_ok=True)
        cmd = ["index_gff", "--index", gff, out]
        if ctx.cfg["annotation"].get("compress_id"):
            cmd.append("--compress-id")
        if run(ctx, cmd):
            sys.exit(f"[ERROR] index_gff failed for {t}")
        print(f"[OK] {t}: {out}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="check BAMs and annotations").set_defaults(fn=step_check)
    p = sub.add_parser("index", help="index_gff for each event type")
    p.add_argument("--event-type", help="only this event type")
    p.set_defaults(fn=step_index)
    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
