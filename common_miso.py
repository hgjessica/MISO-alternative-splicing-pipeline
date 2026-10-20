"""Shared helpers for the MISO step scripts (config, running commands, task selection, SLURM)."""
import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path


class Ctx:
    def __init__(self, cfg, dry_run):
        self.cfg = cfg
        self.dry = dry_run
        self.work = Path(cfg["work_dir"])
        self.threads = str(cfg.get("threads", 4))
        self.miso = cfg.get("miso", {})
        self.types = list(cfg["annotation"]["events"])

    def path(self, *parts):
        return self.work.joinpath(*parts)

    def index_dir(self, etype):
        return self.path("index", etype)

    def miso_dir(self, etype, sample):
        return self.path("miso", etype, sample)

    def compare_dir(self, etype, comparison):
        return self.path("compare", etype, comparison)

    def bf_file(self, etype, comp):
        return self.compare_dir(etype, comp["name"]) / "bayes-factors" / \
            f"{comp['control']}_vs_{comp['treated']}.miso_bf"


def load_config(path):
    with open(path) as fh:
        cfg = json.load(fh)
    for key in ("work_dir", "samples", "annotation"):
        if key not in cfg:
            sys.exit(f"config is missing required key: {key}")
    if not cfg["annotation"].get("events"):
        sys.exit("config: annotation.events must map an event type (SE, A3SS, ...) to a GFF3 file")
    names = [s["name"] for s in cfg["samples"]]
    if len(set(names)) != len(names):
        sys.exit("config: sample names must be unique")
    for c in cfg.get("comparisons", []):
        for k in ("name", "control", "treated"):
            if k not in c:
                sys.exit(f"config: comparison needs '{k}': {c}")
        for k in ("control", "treated"):
            if c[k] not in names:
                sys.exit(f"config: comparison '{c['name']}' refers to unknown sample {c[k]}")
    return cfg


def run(ctx, cmd, env=None, **kw):
    """Run a command (list form, no shell). Returns the exit code."""
    print("$ " + " ".join(shlex.quote(str(c)) for c in cmd), flush=True)
    if ctx.dry:
        return 0
    full_env = dict(os.environ, **(env or {}))
    return subprocess.run([str(c) for c in cmd], env=full_env, **kw).returncode


def need(ctx, *files):
    """Abort with a clear message if an input is missing (skipped in dry run)."""
    for f in files:
        if not ctx.dry and not Path(f).exists():
            sys.exit(f"[ERROR] missing: {f}")


def raise_open_files(limit=10000):
    """MISO opens many files at once; the original jobs ran `ulimit -n 10000`."""
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        want = limit if hard == resource.RLIM_INFINITY else min(limit, hard)
        if soft < want:
            resource.setrlimit(resource.RLIMIT_NOFILE, (want, hard))
    except (ImportError, ValueError, OSError):
        pass


def kept_samples(cfg):
    return [s for s in cfg["samples"] if s.get("keep", True)]


def sample_by_name(cfg, name):
    for s in cfg["samples"]:
        if s["name"] == name:
            return s
    sys.exit(f"unknown sample: {name}")


def read_len(ctx, sample):
    rl = sample.get("read_len", ctx.miso.get("read_len"))
    if rl is None:
        sys.exit(f"no read_len for sample {sample['name']} (set miso.read_len or sample.read_len)")
    return str(rl)


def units(ctx, unit):
    """All tasks of a kind as (object, event_type) pairs, in a fixed order."""
    if unit == "sample":
        return [(s, None) for s in kept_samples(ctx.cfg)]
    if unit == "sample_event":
        return [(s, t) for s in kept_samples(ctx.cfg) for t in ctx.types]
    if unit == "comparison_event":
        return [(c, t) for c in ctx.cfg.get("comparisons", []) for t in ctx.types]
    raise ValueError(unit)


def pick_tasks(ctx, unit, name=None, etype=None):
    """--sample/--comparison NAME and/or --event-type pick tasks locally; under SLURM the array
    task id indexes the task list. Without any selector everything is returned only if the
    caller passes name='ALL'."""
    allu = units(ctx, unit)
    if etype and etype not in ctx.types:
        sys.exit(f"unknown event type {etype}; configured: {', '.join(ctx.types)}")
    if name or etype:
        sel = [(o, t) for o, t in allu
               if (not name or name == "ALL" or o["name"] == name) and (not etype or t == etype)]
        if not sel:
            sys.exit(f"no task matches name={name} event-type={etype}")
        return sel
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    if task is None:
        sys.exit("select a task with --sample/--comparison NAME and/or --event-type "
                 "(or run as a SLURM array task)")
    if int(task) >= len(allu):
        sys.exit(f"array task {task} out of range ({len(allu)} tasks)")
    return [allu[int(task)]]


def task_tmp(ctx, tag):
    """Per-task scratch directory so concurrent array tasks cannot collide."""
    base = Path(ctx.cfg.get("scratch_dir") or ctx.path(".tmp"))
    d = base / f"{tag}_{os.environ.get('SLURM_ARRAY_TASK_ID', 'local')}"
    if not ctx.dry:
        d.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(d)
    return d


def base_parser(doc):
    """Parser with the options every step script shares."""
    ap = argparse.ArgumentParser(description=doc.strip().split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="\n".join(doc.strip().split("\n\n")[1:]))
    ap.add_argument("--config", required=True, help="JSON config (see config.example.json)")
    ap.add_argument("--dry-run", action="store_true", help="print commands, do not run")
    return ap


def add_selectors(p, kind):
    """--sample/--comparison and --event-type options for per-task steps."""
    p.add_argument(f"--{kind}", help=f"{kind} name, or ALL (default: SLURM array task id)")
    p.add_argument("--event-type", help="only this event type (SE, A3SS, ...)")


def add_submit(sub, step_units):
    """Add the `submit` subcommand (SLURM array). step_units: {step name: task unit}."""
    p = sub.add_parser("submit", help="submit a step as a SLURM array")
    p.add_argument("step_name", choices=sorted(step_units))
    p.add_argument("step_args", nargs=argparse.REMAINDER,
                   help="extra args passed to the step")
    return p


def submit_array(ctx, args, script, step_units):
    """Write and submit a SLURM array script (one task per sample/comparison x event type)."""
    sl = ctx.cfg.get("slurm", {})
    n = len(units(ctx, step_units[args.step_name]))
    if n == 0:
        sys.exit("nothing to submit (no samples/comparisons configured)")
    conc = sl.get("array_concurrency", 3)
    extra = " ".join(shlex.quote(a) for a in args.step_args)
    log_dir = ctx.path("slurm_logs")
    text = f"""#!/bin/bash
#SBATCH -p {sl.get('partition', 'PARTITION_NAME')}
#SBATCH -N 1
#SBATCH -n {ctx.threads}
#SBATCH --mem={sl.get('mem', '40G')}
#SBATCH -t {sl.get('time', '144:00:00')}
#SBATCH --export=ALL
#SBATCH -J {args.step_name}
#SBATCH --array=0-{n - 1}%{conc}
#SBATCH -o {log_dir}/{args.step_name}_%A_%a.out
#SBATCH -e {log_dir}/{args.step_name}_%A_%a.err

ulimit -n 10000
export MPLBACKEND=Agg

# the MISO tools live in a python2.7 conda env; the step scripts need a python >= 3.8
source {sl.get('conda_sh', 'Path_to_conda.sh')}
conda activate {sl.get('conda_env', 'miso')}
{sl.get('python', 'python3')} {Path(script).resolve()} --config {Path(args.config).resolve()} {args.step_name} {extra}
"""
    sbatch_file = ctx.path("slurm", f"{args.step_name}.sbatch")
    print(text)
    if ctx.dry:
        return
    log_dir.mkdir(parents=True, exist_ok=True)
    sbatch_file.parent.mkdir(parents=True, exist_ok=True)
    sbatch_file.write_text(text)
    subprocess.run(["sbatch", str(sbatch_file)])
