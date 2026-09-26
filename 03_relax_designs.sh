#!/usr/bin/env bash
# =============================================================================
# 03_relax_designs.sh  —  Amber/GPU relaxation of selected scout designs
# -----------------------------------------------------------------------------
# Re-folds chosen complex(es) and applies ColabFold's --amber --use-gpu-relax,
# producing relaxed_rank_001 models that are directly comparable to your
# Colab-relaxed binders (same minimization path).
#
# Best practice: reuse the cached .a3m from the scout run (MSA_CACHE below).
# Feeding the SAME MSA reproduces the SAME structure, so the only thing that
# changes is the relaxation -> a clean test of whether key residues hold.
# =============================================================================
set -euo pipefail

# ----------------------------- USER SETTINGS ---------------------------------
RF_ROOT="${RF_ROOT:-$HOME/Protein-Design/rf_mpnn_af2}"
COLABFOLD_BIN="$RF_ROOT/tools/localcolabfold/localcolabfold/colabfold-conda/bin"

# complexes.fasta from the scout run / Tier-1 archive (binder:target sequences).
# CONFIRM this path against your actual archive location.
SRC_FASTA="$RF_ROOT/archives/<RUN>/tier1/complexes.fasta"

# Which design jobnames to relax (header names in complexes.fasta). Add more to batch.
DESIGNS=( "design_0_s0" )   # PLACEHOLDER: jobname(s) from complexes.fasta

# OPTIONAL: reuse cached MSAs to skip the server AND reproduce the exact structure.
# Point at the dir holding <jobname>.a3m. Leave "" to fetch MSA fresh (near-identical,
# but not bit-identical to the archived model).
MSA_CACHE=""        # e.g. "$RF_ROOT/workflows/<scout_run>/msas"

# Match these to your ORIGINAL scout fold settings.
MODEL_TYPE="alphafold2_multimer_v3"   # from your filenames
NUM_RECYCLE=3                         # set to whatever the scout used
# Reproduce the exact analyzed model: rank_001 came from model_5, so fold just it.
# (Switch to --num-models 5 / drop --model-order if you'd rather re-rank fresh.)
NUM_MODELS=1
MODEL_ORDER=5
NUM_RELAX=1                           # relax the top-ranked model only
# -----------------------------------------------------------------------------

STAMP=$(date +%Y%m%d_%H%M%S)
RELAX_IN="$RF_ROOT/workflows/relax/${STAMP}/in"
RELAX_OUT="$RF_ROOT/workflows/relax/${STAMP}/out"
LOGDIR="$HOME/Protein-Design/logs"

export PATH="$COLABFOLD_BIN:$PATH"
mkdir -p "$RELAX_IN" "$RELAX_OUT" "$LOGDIR"

# --- Build the relax input: cached a3m if present, else pull from complexes.fasta ---
for d in "${DESIGNS[@]}"; do
  if [ -n "$MSA_CACHE" ] && [ -f "$MSA_CACHE/$d.a3m" ]; then
    cp "$MSA_CACHE/$d.a3m" "$RELAX_IN/$d.a3m"
    echo "[$d] reusing cached MSA (exact-structure relax)"
  else
    # Pull the >jobname record and its (single-line, ':'-separated) complex sequence.
    awk -v id="$d" '$0 ~ "^>"id"([ \t]|$)" {print; getline; print; exit}' \
        "$SRC_FASTA" > "$RELAX_IN/$d.fasta"
    [ -s "$RELAX_IN/$d.fasta" ] || { echo "!! $d not found in $SRC_FASTA"; exit 1; }
    echo "[$d] staged from complexes.fasta (will fetch MSA)"
  fi
done

LOG="$LOGDIR/relax_${STAMP}.log"
echo "Relax -> $RELAX_OUT"
echo "Log    -> $LOG"

nohup bash -c "
  export PATH='$COLABFOLD_BIN:\$PATH'
  colabfold_batch \
    --model-type $MODEL_TYPE \
    --num-recycle $NUM_RECYCLE \
    --num-models $NUM_MODELS --model-order $MODEL_ORDER \
    --amber --use-gpu-relax --num-relax $NUM_RELAX \
    '$RELAX_IN' '$RELAX_OUT'
  echo '=== RELAX DONE ==='
" > "$LOG" 2>&1 &

echo "PID $! detached.  Monitor:  tail -f $LOG   (or your scoutstatus alias)"
