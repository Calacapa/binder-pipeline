# binder-pipeline

A reproducible command-line pipeline for de novo mini-binder design against a
protein target, wrapping RFdiffusion, ProteinMPNN and AlphaFold2-Multimer
(via ColabFold) into three staged drivers with a shared library of stage
functions.

Developed at Calacapa LLC for antivirulence binder design — inhibiting secreted
bacterial virulence factors rather than killing the organism. The code here is
target-agnostic; target structures, designed sequences and campaign results are
not included.

## Design

Two layers, deliberately:

- **`lib/pipeline.sh`** is the single source of truth for how each engine is
  invoked. Every driver calls the same `stage_rfdiffusion`, `stage_mpnn` and
  `stage_af2` functions, so a change to how a tool is called happens in exactly
  one place.
- **The numbered drivers** decide only *which* stages run and *with what
  numbers*. They contain no engine-specific arguments.

This split is what makes runs comparable across campaigns: an explore run and a
refine run differ in sampling depth, not in how the models were called.

Every run writes to a timestamped folder under `$RF_ROOT/runs/`, so nothing is
ever overwritten and any result can be traced back to the parameters that
produced it (`run.info` in each run folder).

## The three stages

### `01_explore_new_epitope.sh` — scout a new epitope

Full funnel: RFdiffusion backbones → ProteinMPNN sequences → AF2-Multimer
co-folding → metric filter. Two presets:

| mode | backbones | seqs/backbone | AF2 jobs | use |
|---|---|---|---|---|
| `light` | 20 | 4 | ~80 | scout a new patch cheaply |
| `standard` | 100 | 8 | ~800 | commit to a patch |

```bash
./01_explore_new_epitope.sh --mode light
./01_explore_new_epitope.sh --mode standard --hotspots A30,A33,A34
./01_explore_new_epitope.sh --mode light --backbones 12
```

### `02_refine_top_binder.sh` — deep-dive on backbones you like

Skips RFdiffusion entirely and re-runs the *same* MPNN and AF2 stages on
backbones you already have, with many more sequences each (low volume, high
depth) and a stricter default gate.

```bash
./02_refine_top_binder.sh --in ./my_top_hits --seqs 48 --iptm 0.85
```

Input PDBs must follow the convention explore produces: chain A = binder,
chain B = target.

### `03_relax_designs.sh` — Amber/GPU relaxation

Re-folds selected complexes with ColabFold's `--amber --use-gpu-relax`,
producing relaxed models comparable to the originals. Reusing the cached `.a3m`
from the original run reproduces the same structure, so relaxation is the only
variable — a clean test of whether key interface residues hold.

## Filtering

Designs are gated on three AF2-Multimer metrics, at field-standard defaults:

| metric | default | meaning |
|---|---|---|
| ipTM | ≥ 0.80 | predicted interface accuracy |
| interface PAE | ≤ 10.0 Å | predicted error across the binder–target interface |
| binder pLDDT | ≥ 80.0 | per-residue confidence in the binder chain |

`02` tightens these to 0.85 / 8.0 / 82.0 by default. All are overridable per run
via `--iptm`, `--pae`, `--plddt`.

`lib/af2_glue.py` does the two jobs that sit between engines: `make-fasta`
assembles binder:target complex FASTAs from ProteinMPNN output with a manifest,
and `filter` parses AF2 output, computes interface PAE and binder pLDDT, writes
a ranked `metrics.csv` and copies passing structures to `hits/`.

## Requirements

Not a turnkey install. Assumes a GPU host with:

- RFdiffusion, in a conda env (`$RF_ROOT/envs/SE3nv`)
- ProteinMPNN, in a conda env (`$RF_ROOT/envs/mpnn`)
- LocalColabFold providing `colabfold_batch` on `PATH`
- conda available to non-interactive shells

Set `RF_ROOT` to your install root (defaults to
`$HOME/Protein-Design/rf_mpnn_af2`). Drivers locate `lib/` relative to
themselves, so the repo can live anywhere.

Developed and run on remote GPU instances driven from JupyterLab.

## Placeholders

The `CONFIG` blocks ship with placeholder values — target range, hotspot
residues, design jobnames. Replace them with your own; they are marked
`PLACEHOLDER` in the source.

## Citation

See `CITATION.cff`.

## License

Apache License 2.0.
