#!/usr/bin/env bash
# ======================================================================
# 01_explore_new_epitope.sh
#   Full funnel for a BRAND-NEW epitope: RFdiffusion -> MPNN -> AF2 -> filter.
#   Two presets via --mode:
#     light     ~20 backbones x 4 seqs  -> 80 AF2   (scout a new patch)
#     standard  ~100 backbones x 8 seqs -> 800 AF2  (commit to the patch)
#
# Usage:
#   ./01_explore_new_epitope.sh --mode light
#   ./01_explore_new_epitope.sh --mode standard --hotspots A30,A33,A34
#   ./01_explore_new_epitope.sh --mode light --backbones 12   # override a number
#
# Edit the CONFIG block once per target. Everything else is shared lib.
# ======================================================================
set -euo pipefail

# ---------------- CONFIG (edit per target) ----------------------------
RF_ROOT="${RF_ROOT:-$HOME/Protein-Design/rf_mpnn_af2}"
TARGET_PDB="$RF_ROOT/shared/targets/target.pdb"
TARGET_CHAIN="A"          # chain of the TARGET in TARGET_PDB
TARGET_RANGE="1-150"      # PLACEHOLDER: residues of the target to present
BINDER_LEN="70-100"       # binder length range RFdiffusion samples
HOTSPOTS="A30,A33,A34"    # PLACEHOLDER: epitope patch (target-chain numbering)

# AF2 gate (field-standard). Tune here or via flags.
IPTM_MIN="0.80"
PAE_MAX="10.0"
PLDDT_MIN="80.0"
# ----------------------------------------------------------------------

MODE="light"
BACKBONES=""; SEQS=""     # empty = use the mode preset

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)      MODE="$2"; shift 2;;
    --backbones) BACKBONES="$2"; shift 2;;
    --seqs)      SEQS="$2"; shift 2;;
    --hotspots)  HOTSPOTS="$2"; shift 2;;
    --target)    TARGET_PDB="$2"; shift 2;;
    --iptm)      IPTM_MIN="$2"; shift 2;;
    --pae)       PAE_MAX="$2"; shift 2;;
    --plddt)     PLDDT_MIN="$2"; shift 2;;
    *) echo "unknown arg: $1"; exit 1;;
  esac
done

# preset -> numbers (overridable by flags above)
case "$MODE" in
  light)    : "${BACKBONES:=20}";  : "${SEQS:=4}";;
  standard) : "${BACKBONES:=100}"; : "${SEQS:=8}";;
  *) echo "mode must be 'light' or 'standard'"; exit 1;;
esac

source "$(dirname "$(readlink -f "$0")")/lib/pipeline.sh"

# Timestamped run folder so nothing is ever overwritten.
RUN="$RF_ROOT/runs/explore_${MODE}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN"
CONTIG="${TARGET_CHAIN}${TARGET_RANGE}/0 ${BINDER_LEN}"

log "EXPLORE [$MODE]  backbones=$BACKBONES  seqs/bb=$SEQS  -> $RUN"
echo "target=$TARGET_PDB | contig=[$CONTIG] | hotspots=$HOTSPOTS" | tee "$RUN/run.info"

stage_rfdiffusion "$TARGET_PDB" "$CONTIG" "$HOTSPOTS" "$BACKBONES" "$RUN/rfd"
stage_mpnn        "$RUN/rfd" "$SEQS" "$RUN/mpnn" "A"
stage_af2         "$RUN/mpnn/seqs" "$RUN/af2" "$IPTM_MIN" "$PAE_MAX" "$PLDDT_MIN"

log "DONE. Ranked hits: $RUN/af2/metrics.csv"
column -s, -t "$RUN/af2/metrics.csv" | head -20 || cat "$RUN/af2/metrics.csv"
