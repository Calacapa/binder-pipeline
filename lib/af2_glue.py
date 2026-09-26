#!/usr/bin/env python3
"""
af2_glue.py  --  the two glue steps between MPNN and AF2, and after AF2.

Two subcommands, no third-party deps (stdlib + numpy only):

  make-fasta   MPNN seqs/*.fa  ->  one combined ColabFold complex FASTA
               (binder:target per design) + a manifest of binder lengths.

  filter       ColabFold output dir + manifest  ->  ranked metrics.csv
               + copies of passing complex PDBs into a hits/ folder.

Design conventions (match RFdiffusion binder output + ProteinMPNN):
  * RFdiffusion writes chain A = designed binder, chain B = target.
  * ProteinMPNN concatenates chains in the .fa with '/', in order A/B,
    so chunk[0] = binder, chunk[-1] = target. The TARGET sequence is read
    from the *native* (first) record, which is the fixed, true target.
  * ColabFold represents a complex as SEQ1:SEQ2 in one FASTA entry.
"""
import argparse, csv, glob, json, os, re, shutil, sys

# numpy is part of every AF2/colabfold env; import lazily so make-fasta works
# even in a bare python.
def _np():
    import numpy as np
    return np


# ----------------------------------------------------------------------
# make-fasta : MPNN .fa  ->  ColabFold complex FASTA + manifest
# ----------------------------------------------------------------------
def read_fasta(path):
    """Return list of (header, sequence) preserving order."""
    recs, hdr, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr, seq = line[1:], []
            else:
                seq.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    return recs


def make_fasta(args):
    fa_files = sorted(glob.glob(os.path.join(args.mpnn_seqs, "*.fa")) +
                      glob.glob(os.path.join(args.mpnn_seqs, "*.fasta")))
    if not fa_files:
        sys.exit(f"[make-fasta] no .fa files in {args.mpnn_seqs}")

    out_fa = open(args.out_fasta, "w")
    manifest = open(args.manifest, "w", newline="")
    mw = csv.writer(manifest)
    mw.writerow(["name", "binder_len", "target_len", "source_fa"])

    n_written = 0
    for fa in fa_files:
        stem = re.sub(r"\.(fa|fasta)$", "", os.path.basename(fa))
        recs = read_fasta(fa)
        if not recs:
            continue
        # native (first) record holds the true, fixed target sequence
        native_chunks = recs[0][1].split("/")
        if len(native_chunks) < 2:
            print(f"[make-fasta] WARN {stem}: single-chain .fa, skipping "
                  f"(expected binder/target).", file=sys.stderr)
            continue
        target_seq = native_chunks[args.target_index].strip()

        # designed records (skip native at index 0)
        for i, (hdr, seq) in enumerate(recs[1:], start=1):
            chunks = seq.split("/")
            binder_seq = chunks[args.binder_index].strip()
            if not binder_seq:
                continue
            name = f"{stem}_s{i}"
            out_fa.write(f">{name}\n{binder_seq}:{target_seq}\n")
            mw.writerow([name, len(binder_seq), len(target_seq),
                         os.path.basename(fa)])
            n_written += 1

    out_fa.close()
    manifest.close()
    print(f"[make-fasta] wrote {n_written} complex entries -> {args.out_fasta}")
    if n_written == 0:
        sys.exit("[make-fasta] nothing written; check MPNN output.")


# ----------------------------------------------------------------------
# filter : ColabFold output  ->  ranked metrics.csv + hits/
# ----------------------------------------------------------------------
def load_manifest(path):
    m = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            m[row["name"]] = int(row["binder_len"])
    return m


def find_top_score_json(out_dir, name):
    """ColabFold writes <name>_scores_rank_001_*.json (rank 1 = best)."""
    hits = sorted(glob.glob(os.path.join(out_dir, f"{name}_scores_rank_001*.json")))
    if not hits:
        hits = sorted(glob.glob(os.path.join(out_dir, f"{name}_scores_rank_*001*.json")))
    return hits[0] if hits else None


def find_top_pdb(out_dir, name):
    hits = sorted(glob.glob(os.path.join(out_dir, f"{name}_*rank_001*.pdb")))
    return hits[0] if hits else None


def interface_pae(pae, binder_len):
    """Mean PAE of the two off-diagonal (binder<->target) blocks."""
    np = _np()
    pae = np.asarray(pae, dtype=float)
    Lb = binder_len
    bt = pae[:Lb, Lb:]
    tb = pae[Lb:, :Lb]
    vals = []
    if bt.size:
        vals.append(bt.mean())
    if tb.size:
        vals.append(tb.mean())
    return float(np.mean(vals)) if vals else float("nan")


def filter_results(args):
    np = _np()
    manifest = load_manifest(args.manifest)
    os.makedirs(args.hits_dir, exist_ok=True)
    rows = []

    for name, Lb in manifest.items():
        sj = find_top_score_json(args.af2_dir, name)
        if not sj:
            print(f"[filter] no scores json for {name}", file=sys.stderr)
            continue
        with open(sj) as fh:
            d = json.load(fh)

        plddt = np.asarray(d.get("plddt", []), dtype=float)
        pae = d.get("pae", d.get("predicted_aligned_error"))
        iptm = float(d.get("iptm", d.get("ptm", float("nan"))))

        binder_plddt = float(plddt[:Lb].mean()) if plddt.size >= Lb else float("nan")
        pae_int = interface_pae(pae, Lb) if pae is not None else float("nan")

        passed = (iptm >= args.iptm_min and
                  pae_int <= args.pae_max and
                  binder_plddt >= args.plddt_min)

        rows.append({
            "name": name, "binder_len": Lb,
            "iptm": round(iptm, 4),
            "pae_interaction": round(pae_int, 3),
            "binder_plddt": round(binder_plddt, 2),
            "pass": int(passed),
        })

        if passed:
            pdb = find_top_pdb(args.af2_dir, name)
            if pdb:
                shutil.copy(pdb, os.path.join(args.hits_dir, os.path.basename(pdb)))

    # rank: passing first, then by iptm desc, then pae asc
    rows.sort(key=lambda r: (-r["pass"], -r["iptm"], r["pae_interaction"]))
    with open(args.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "binder_len", "iptm",
                                           "pae_interaction", "binder_plddt", "pass"])
        w.writeheader()
        w.writerows(rows)

    n_pass = sum(r["pass"] for r in rows)
    print(f"[filter] scored {len(rows)} designs, {n_pass} passed "
          f"(iptm>={args.iptm_min}, pae_int<={args.pae_max}, plddt>={args.plddt_min})")
    print(f"[filter] metrics -> {args.out_csv} ; hits -> {args.hits_dir}/")


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    mf = sub.add_parser("make-fasta")
    mf.add_argument("--mpnn-seqs", required=True, help="MPNN out_folder/seqs dir")
    mf.add_argument("--out-fasta", required=True)
    mf.add_argument("--manifest", required=True)
    mf.add_argument("--binder-index", type=int, default=0,
                    help="chunk index of binder in MPNN '/' split (default 0 = chain A)")
    mf.add_argument("--target-index", type=int, default=-1,
                    help="chunk index of target (default -1 = last chain)")
    mf.set_defaults(func=make_fasta)

    fl = sub.add_parser("filter")
    fl.add_argument("--af2-dir", required=True, help="colabfold output dir")
    fl.add_argument("--manifest", required=True)
    fl.add_argument("--out-csv", required=True)
    fl.add_argument("--hits-dir", required=True)
    fl.add_argument("--iptm-min", type=float, default=0.80)
    fl.add_argument("--pae-max", type=float, default=10.0)
    fl.add_argument("--plddt-min", type=float, default=80.0)
    fl.set_defaults(func=filter_results)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
