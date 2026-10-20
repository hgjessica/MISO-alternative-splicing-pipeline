# MISO-alternative-splicing-pipeline
Python scripts that run [MISO](http://genes.mit.edu/burgelab/miso/) (Mixture of Isoforms, `misopy` 0.5.4) on aligned short-read RNA-seq: estimate the percent spliced in (Psi) of annotated alternative-splicing events per sample, compare two conditions with Bayes factors, and keep the significant events. Runs locally or as SLURM job arrays.

| Script | Covers |
|---|---|
| `01_prepare_events.py` | `check` inputs, `index_gff` for each event type (SE, A3SS, A5SS, MXE, RI) |
| `02_run_miso.py` | `miso --run` per sample and event type (+ optional paired-end `insert-len`) |
| `03_summarize_compare_filter.py` | `summarize_miso`, `compare_miso`, delta-Psi / Bayes factor filter |
| `04_sashimi_plot.py` | optional `sashimi_plot` for chosen events, Bayes factor distribution |

`common.py` holds the shared helpers and must stay in the same folder. Paths, sample names and
cluster settings live in one JSON config, so the code contains nothing project-specific.

```
BAM (STAR) + event GFF3 ──index_gff──▶ index/<TYPE>
index + BAM ──miso --run──▶ miso/<TYPE>/<sample> ──compare_miso──▶ compare/<TYPE>/<comparison>/*.miso_bf
                                                      └─ filter ─▶ filtered/<TYPE>/<comparison>.significant.tsv
```

The STAR (or other splice-aware) alignment is not part of this pipeline; point `samples[].bam`
at sorted, indexed BAMs.

## 1. Requirements

- **Two Pythons.** MISO is a Python 2.7 program (the original analysis ran in python 2.7.18,
  `misopy 0.5.4`, numpy 1.16, scipy 1.2.3, matplotlib 2.2.5, pysam 0.20.0, samtools 1.5,
  bedtools 2.26.0). The step scripts here are Python >= 3.8, standard library only, and just call
  the MISO commands. So: activate the `miso` conda env (puts `miso`, `index_gff`,
  `compare_miso`, `summarize_miso`, `sashimi_plot`, `samtools` on `PATH`) and run the scripts
  with a separate `python3`. That env's own python is 2.7, so set `slurm.python` to a python3
  (an absolute path, or one from `module load`) if `python3` is not on `PATH` after activation.
- Optionally SLURM (`sbatch`) for the `submit` command.

```bash
conda env create -f environment.yml     # python 2.7 MISO env; see notes inside the file
conda activate miso
```

`environment.yml` records the versions found in the env that was actually used. It has not
been re-created from this file, and the original env was built from a MISO source checkout, so
if `pip install misopy==0.5.4` fails, install MISO from its GitHub source.

You also need:
- event annotations as MISO GFF3 files, one per event type (`SE`, `A3SS`, `A5SS`, `MXE`, `RI`),
  with `chr`-prefixed chromosome names. See **Event annotation files** below for how to obtain them.
- sorted, indexed BAM files with headers whose chromosome names match the GFF3 (`chr1` vs `1`).

## 2. Configure

Copy `config.example.json` to `config.json` and fill in every `Path_to_...` value.

| Key | Meaning |
|---|---|
| `work_dir` | All outputs go here |
| `threads` | CPUs per SLURM task (`-n`) |
| `scratch_dir` | Optional. Per-task temp dir |
| `slurm` | Partition, memory, time, array concurrency, conda setup, `python` (python3 used for the scripts). Defaults 40G / 144 h as in the original jobs |
| `annotation.events` | `{event type: GFF3}`. Every type listed is indexed and run; remove types you do not need. `compress_id` adds `--compress-id` to `index_gff` |
| `miso` | `read_len` (default for samples), `paired_end` (`{"mean":..,"sd":..}` or null), `overhang_len`, `settings_filename`, `num_processes` (`-p`), `no_filter_events`, `prefilter`, `pass_event_type` |
| `samples` | `name`, `bam`, optional `read_len` (overrides `miso.read_len`), optional `paired_end`, `keep` |
| `comparisons` | `name`, `control`, `treated` (sample names). Sample 1 is the control, so `diff = PSI_control - PSI_treated` |
| `filter` | `delta_psi`, `bayes_factor`, `num_inc`, `num_exc`, `num_sum_inc_exc`, `num_total`, `apply_both` |
| `plots` | `top_n` (plot the N best significant events), `events` (`{"type","event"}` list), `settings` (overrides of the plot parameters) |

Parameters carried over from the original runs: `miso --run <index> <bam> --output-dir <dir>
--read-len N` with a per-dataset read length (75, 100 and 150 were used), no other flags (so
MISO's defaults: single-end mode, overhang 1, 4 processes, event filter of at least 20 reads,
settings from MISO's own `miso_settings.txt`: burn-in 500, lag 10, 5000 iterations, 6 chains),
`compare_miso --compare-samples <control> <treated> <outdir>`, and the plot settings
(`intron_scale 30`, `exon_scale 4`, `font_size 6`, `resolution .5`, `posterior_bins 40`,
`gene_posterior_ratio 5`, `bf_thresholds [0,1,2,5,10,20]`).

**The filter defaults are not from the original analysis.** The original scripts never called
`filter_events`; events were classified in notebooks from the sign of `diff` and the Bayes
factor. The defaults (`delta_psi 0.20`, `bayes_factor 10`, `num_inc 1`, `num_exc 1`,
`num_sum_inc_exc 10`) are the example values from the MISO manual. Set your own.

Keep `config.json` out of version control if it contains real names or paths.

## 3. Run

Add `--dry-run` before the subcommand to print the commands without running anything.
Per-task steps take `--sample`/`--comparison NAME` (or `ALL`) and `--event-type`; under SLURM
they use the array task id instead (one task per sample x event type, or per comparison x event
type).

```bash
P1="python3 01_prepare_events.py --config config.json"
P2="python3 02_run_miso.py --config config.json"
P3="python3 03_summarize_compare_filter.py --config config.json"
P4="python3 04_sashimi_plot.py --config config.json"
```

| Step | Command | Output |
|---|---|---|
| check inputs | `$P1 check` | report only |
| index events (once) | `$P1 index` | `index/<TYPE>/` |
| Psi per sample | `$P2 submit run` (or `$P2 run --sample S --event-type SE`) | `miso/<TYPE>/<sample>/<chrom>/*.miso` |
| insert size (optional, paired-end) | `$P2 insert-len --sample S` | `insert_len/<sample>/` |
| summaries | `$P3 submit summarize` | `summary/<TYPE>/<sample>/summary/<sample>.miso_summary` |
| compare conditions | `$P3 submit compare` | `compare/<TYPE>/<comparison>/bayes-factors/<control>_vs_<treated>.miso_bf` |
| filter | `$P3 filter --comparison ALL` | `filtered/<TYPE>/<comparison>.{all_events,significant}.tsv`, `filtered/filter_summary.tsv` |
| sashimi plots (optional) | `$P4 submit plot` | `sashimi/<TYPE>/<comparison>/<event>.pdf` |
| Bayes factor distribution | `$P4 bf-dist --comparison ALL` | same folder |

Notes:
- `submit STEP` writes `slurm/STEP.sbatch`, prints it, and calls `sbatch`. The sbatch script
  runs `ulimit -n 10000`, activates the conda env and starts the script with `slurm.python`.
  Test one task locally first.
- `compare` needs the `run` outputs of both samples; wait for the previous array to finish.
- The significant table has one row per event with `psi_control`, `psi_treated`, `diff`,
  `bayes_factor`, `direction` (`control_higher` / `treated_higher`) and the inclusion /
  exclusion / both-isoform read counts of each sample. Psi is the fraction of the first isoform
  (MISO's convention: inclusion isoform) so a positive `diff` means more inclusion in the control.
- The filter mirrors `filter_events`: an event passes if `|diff|` and the Bayes factor reach the
  thresholds and the read-count rules hold in either sample (both with `apply_both`). Events
  with more than two isoforms are skipped and counted in `filter_summary.tsv`.
- Sashimi coverages (mapped reads, used to scale read density) are computed with
  `samtools view -c -F 260`; set `samples[].mapped_reads` to use another definition.

## Event annotation files

`annotation.events` in the config points to one GFF3 per event type, for example
`Path_to_annotation_dir/SE.hg38.gff3`. These files only describe public gene annotation, so they
can be rebuilt or downloaded; they are not project data. Each file has `gene` / `mRNA` (`.A`, `.B`)
records per event, with IDs of the form `chr1:start:end:strand@chr1:...`, and `chr`-prefixed names.

Ways to get them (check each source yourself before relying on it; none was run for this README):

1. **Prebuilt MISO annotations.** The MISO documentation
   ([annotation page](https://miso.readthedocs.io/en/dev/annotation.html)) lists downloadable event
   sets for hg18/hg19 and mm9/mm10. They are for hg19, not hg38, and were themselves lifted over from
   older assemblies. The MISO users list reports that there is no official hg38 set
   ([thread](https://mailman.mit.edu/pipermail/miso-users/Week-of-Mon-20210830/000419.html)), and a
   link in the docs has been reported as broken.
2. **Lift over hg19 to hg38.** Convert the GFF3 coordinates with UCSC `liftOver`, then re-check that
   exon boundaries still make sense and that no event lost a part.
3. **Build from an hg38 gene annotation.** Create events from a UCSC gene table (for example
   `ensGene`) or GENCODE/Ensembl GTF, then format them as MISO GFF3. The original analysis did this
   with an in-house, unpublished tool using the `commonshortest` flanking rule (the most common
   shortest form of the flanking exons), so there is no public download of the exact files used here.
   Keep the annotation release, genome build and event-building rule in your methods.
4. **Bioconductor `alternativeSplicingEvents.hg38`.** Provides hg38 event definitions compiled from
   MISO, VAST-TOOLS, SUPPA and rMATS
   ([manual](https://bioconductor.statistik.tu-dortmund.de/packages/3.21/data/annotation/manuals/alternativeSplicingEvents.hg38/man/alternativeSplicingEvents.hg38.pdf)).
   It is an R data package, so you would have to convert its events to MISO GFF3 yourself; I have not
   verified a direct export.

Whichever you use:
- Use the **same genome build and chromosome naming** as the BAMs (`01_prepare_events.py check` tests this).
- Prefer restricting the annotation to primary chromosomes before `index`. Alt/fix contigs create
  hundreds of extra index folders and slow `miso --run` start-up.
- Record the annotation source and version next to your results.
