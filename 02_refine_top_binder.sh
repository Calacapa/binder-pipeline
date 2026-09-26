#!/usr/bin/env bash
# ======================================================================
# 02_refine_top_binder.sh
#   DEEP-DIVE on 1-2 backbones you already like (e.g. hits from explore).
#   Skips RFdiffusion entirely; reuses the SAME MPNN + AF2 + filter stages,
#   just with many more sequences per backbone (low volume, high depth).
#
# Usage:
#   ./02_refine_top_binder.sh --in <dir_of_pdbs>
#   ./02_refine_top_binder.sh --in ./my_top_hits --seqs 48 --iptm 0.85
#
# <dir_of_pdbs> must contain complex backbones with chain A = binder,
# chain B = target (i.e. the same convention explore produces).
# ======================================================================
set -euo pipefail

# ---------------- CONFIG ----------------------------------------------
IN_DIR=""                 # folder of backbone PDBs to refine (required)
SEQS="48"                 # many sequence variants per backbone
# AF2 gate slightly stricter by default for a refine pass:
IPTM_MIN="0.85"
PAE_MAX="8.0"
PLDDT_MIN="82.0"
# ----------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in)     IN_DIR="$2"; shift 2;;
    --seqs)   SEQS="$2"; shift 2;;
    --iptm)   IPTM_MIN="$2"; shift 2;;
    --pae)    PAE_MAX="$2"; shift 2;;
    --plddt)  PLDDT_MIN="$2"; shift 2;;
    *) echo "unknown arg: $1"; exit 1;;
  esac
done

[[ -z "$IN_DIR" ]] && { echo "need --in <dir_of_backbone_pdbs>"; exit 1; }
[[ -d "$IN_DIR" ]] || { echo "not a directory: $IN_DIR"; exit 1; }

source "$(dirname "$(readlink -f "$0")")/lib/pipeline.sh"

RUN="$RF_ROOT/runs/refine_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN"
cp "$IN_DIR"/*.pdb "$RUN/" 2>/dev/null || true
n_pdb=$(ls "$RUN"/*.pdb 2>/dev/null | wc -l)
[[ "$n_pdb" -eq 0 ]] && { echo "no .pdb files in $IN_DIR"; exit 1; }

log "REFINE  backbones=$n_pdb  seqs/bb=$SEQS  -> $RUN"
echo "in=$IN_DIR | seqs=$SEQS | gate iptm>=$IPTM_MIN pae<=$PAE_MAX plddt>=$PLDDT_MIN" \
  | tee "$RUN/run.info"

# No RFdiffusion stage -- feed the provided backbones straight to MPNN.
stage_mpnn "$RUN" "$SEQS" "$RUN/mpnn" "A"
stage_af2  "$RUN/mpnn/seqs" "$RUN/af2" "$IPTM_MIN" "$PAE_MAX" "$PLDDT_MIN"

log "DONE. Ranked refined variants: $RUN/af2/metrics.csv"
column -s, -t "$RUN/af2/metrics.csv" | head -30 || cat "$RUN/af2/metrics.csv"
