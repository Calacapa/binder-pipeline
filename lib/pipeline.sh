#!/usr/bin/env bash
# ======================================================================
# lib/pipeline.sh  --  shared stage functions for the binder pipeline.
# Sourced by 01_explore_new_epitope.sh and 02_refine_top_binder.sh.
# This is the SINGLE SOURCE OF TRUTH for how each engine is called;
# the drivers only decide which stages run and with what numbers.
# ======================================================================

# Resolve roots (respects activate.sh if already sourced).
: "${RF_ROOT:=$HOME/Protein-Design/rf_mpnn_af2}"
SE3_ENV="$RF_ROOT/envs/SE3nv"          # RFdiffusion
MPNN_ENV="$RF_ROOT/envs/mpnn"          # ProteinMPNN (also has numpy for glue)
RFD="$RF_ROOT/tools/RFdiffusion/scripts/run_inference.py"
MPNN_DIR="$RF_ROOT/tools/ProteinMPNN"
GLUE="${GLUE:-$(dirname "${BASH_SOURCE[0]}")/af2_glue.py}"

# Make conda available even in a non-interactive script.
source "$RF_ROOT/miniconda3/etc/profile.d/conda.sh"

log() { echo -e "\n[$(date +%H:%M:%S)] === $* ==="; }

# ----------------------------------------------------------------------
# stage_rfdiffusion  <target_pdb> <contig> <hotspots> <num_designs> <out_dir>
#   Generates binder backbones. Output: <out_dir>/design_*.pdb
#   NOTE the HYDRA gotcha from BUILD_LOG: output_prefix MUST be absolute.
# ----------------------------------------------------------------------
stage_rfdiffusion() {
  local target_pdb="$1" contig="$2" hotspots="$3" ndes="$4" out_dir="$5"
  mkdir -p "$out_dir"
  local prefix; prefix="$(readlink -f "$out_dir")/design"   # absolute
  log "RFdiffusion: $ndes backbones | contig=$contig | hotspots=$hotspots"
  conda activate "$SE3_ENV"
  python "$RFD" \
    inference.input_pdb="$(readlink -f "$target_pdb")" \
    inference.output_prefix="$prefix" \
    "contigmap.contigs=[$contig]" \
    "ppi.hotspot_res=[$hotspots]" \
    inference.num_designs="$ndes" \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
  conda deactivate
  echo "[rfdiffusion] backbones -> $out_dir"
}

# ----------------------------------------------------------------------
# stage_mpnn  <pdb_dir> <seqs_per_backbone> <out_dir> [design_chains]
#   Designs sequences for the binder chain, target held fixed as context.
#   Output: <out_dir>/seqs/*.fa  (chains joined by '/', binder first)
# ----------------------------------------------------------------------
stage_mpnn() {
  local pdb_dir="$1" nseq="$2" out_dir="$3" design_chains="${4:-A}"
  mkdir -p "$out_dir"
  log "ProteinMPNN: $nseq seq/backbone | design chain(s)=$design_chains | T=0.1"
  conda activate "$MPNN_ENV"
  python "$MPNN_DIR/helper_scripts/parse_multiple_chains.py" \
    --input_path="$pdb_dir" --output_path="$out_dir/parsed.jsonl"
  python "$MPNN_DIR/helper_scripts/assign_fixed_chains.py" \
    --input_path="$out_dir/parsed.jsonl" \
    --output_path="$out_dir/assigned.jsonl" \
    --chain_list "$design_chains"
  python "$MPNN_DIR/protein_mpnn_run.py" \
    --jsonl_path "$out_dir/parsed.jsonl" \
    --chain_id_jsonl "$out_dir/assigned.jsonl" \
    --out_folder "$out_dir" \
    --num_seq_per_target "$nseq" \
    --sampling_temp "0.1" \
    --batch_size 1
  conda deactivate
  echo "[mpnn] sequences -> $out_dir/seqs"
}

# ----------------------------------------------------------------------
# stage_af2  <mpnn_seqs_dir> <out_dir> <iptm_min> <pae_max> <plddt_min>
#   Builds complex FASTAs, co-folds binder+target, filters to hits.
#   Output: <out_dir>/metrics.csv (ranked) + <out_dir>/hits/*.pdb
# ----------------------------------------------------------------------
stage_af2() {
  local seqs_dir="$1" out_dir="$2" iptm="$3" pae="$4" plddt="$5"
  mkdir -p "$out_dir"
  conda activate "$MPNN_ENV"   # has numpy for the glue
  log "AF2 prep: building complex FASTAs from MPNN output"
  python "$GLUE" make-fasta \
    --mpnn-seqs "$seqs_dir" \
    --out-fasta "$out_dir/complexes.fasta" \
    --manifest "$out_dir/manifest.csv"
  conda deactivate

  log "AF2: co-folding complexes with colabfold_batch (GPU)"
  colabfold_batch "$out_dir/complexes.fasta" "$out_dir/af2_out"

  conda activate "$MPNN_ENV"
  log "Filter: iptm>=$iptm  pae_int<=$pae  binder_plddt>=$plddt"
  python "$GLUE" filter \
    --af2-dir "$out_dir/af2_out" \
    --manifest "$out_dir/manifest.csv" \
    --out-csv "$out_dir/metrics.csv" \
    --hits-dir "$out_dir/hits" \
    --iptm-min "$iptm" --pae-max "$pae" --plddt-min "$plddt"
  conda deactivate
  echo "[af2] results -> $out_dir/metrics.csv  (hits in $out_dir/hits/)"
}
