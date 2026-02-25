#!/usr/bin/env python3
"""AlphaScan — AlphaFold3 result analysis from zip archives."""

import argparse
import fnmatch
import io
import json
import re
import sys
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Wedge, PathPatch
from matplotlib.path import Path as MplPath

from Bio.Data import IUPACData
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.PDBExceptions import PDBConstructionWarning

PLOTS_OUTPUT_DIR = Path("final_interaction_plots")
CORE_SUMMARY_FILENAME = "afm_core_scores.csv"

RECEPTOR_CHAIN = "A"
LIGAND_CHAIN = "B"
CONTACT_DISTANCE_CUTOFF = 4.5
SUMMARY_PAE_THRESHOLD = 5.0

parser_bio = MMCIFParser(QUIET=True)
AA_3_TO_1 = {k.upper(): v for k, v in IUPACData.protein_letters_3to1.items()}


# ── Zip access layer ─────────────────────────────────────────────────

class ZipJobSource:
    """Wraps a ZipFile to read AF3 results by job name."""

    _NON_JOB_DIRS = {"msas", "templates", "__MACOSX"}

    def __init__(self, zip_path: Path):
        self.zip_path = zip_path
        self.zf = zipfile.ZipFile(zip_path, "r")
        self._index: dict[str, list[str]] = {}
        root_files: list[str] = []
        for name in self.zf.namelist():
            parts = name.split("/")
            if len(parts) >= 2 and parts[0]:
                self._index.setdefault(parts[0], []).append(name)
            elif len(parts) == 1 and parts[0]:
                root_files.append(name)

        real_jobs = [k for k in self._index if k not in self._NON_JOB_DIRS]
        if not real_jobs and root_files:
            job_name = None
            for f in root_files:
                m = re.match(r"^(.+?)_model_\d+\.cif$", f)
                if m:
                    job_name = m.group(1)
                    break
            if job_name is None:
                job_name = re.sub(r"\s*\(\d+\)$", "", zip_path.stem)
            all_members = root_files[:]
            for members in self._index.values():
                all_members.extend(members)
            self._index = {job_name: all_members}

    def job_names(self) -> list[str]:
        return sorted(self._index.keys())

    def find(self, job_name: str, pattern: str) -> str | None:
        for member in self._index.get(job_name, []):
            basename = member.split("/")[-1]
            if fnmatch.fnmatch(basename, pattern):
                return member
        return None

    def read_text(self, member_name: str) -> str:
        return self.zf.read(member_name).decode("utf-8")

    def close(self):
        self.zf.close()


# ── Core helpers ─────────────────────────────────────────────────────

def load_pae(src: ZipJobSource, member_name: str):
    d = json.loads(src.read_text(member_name))
    m = d.get("pae") or d.get("predicted_aligned_error")
    if m is None:
        raise KeyError(f"No PAE matrix in {member_name}")
    return np.array(m, dtype=float), d


def extract_iptm(d: dict):
    for k in ("iptm", "ipTM", "iptm_score", "iptmScore"):
        if k in d:
            try:
                return float(d[k])
            except Exception:
                pass
    return np.nan


def read_model0_iptm_fast(src: ZipJobSource, job_name: str) -> float:
    conf = src.find(job_name, "*_summary_confidences_0.json")
    if conf:
        try:
            return float(json.loads(src.read_text(conf)).get("iptm", np.nan))
        except Exception:
            pass
    full = src.find(job_name, "*_full_data_0.json")
    if full:
        try:
            d = json.loads(src.read_text(full))
            val = extract_iptm(d)
            return float(val) if not np.isnan(val) else np.nan
        except Exception:
            pass
    return np.nan


def residue_index_map(model):
    idx_map, idx = {}, 0
    for ch in sorted(model, key=lambda c: c.id):
        for r in sorted(ch.get_residues(), key=lambda r: r.id[1]):
            if r.id[0] == " ":
                idx_map[(ch.id, r.id)] = idx
                idx += 1
    return idx_map


def get_sequences_and_plddt(structure):
    data = {"receptor_seq": "", "ligand_seq": "", "receptor_plddt": [], "ligand_plddt": []}
    for ch in sorted(structure[0], key=lambda c: c.id):
        seq, plddt = [], []
        residues = sorted([r for r in ch if r.id[0] == " "], key=lambda r: r.id[1])
        for r in residues:
            seq.append(AA_3_TO_1.get(r.get_resname().upper(), "X"))
            plddt.append(r["CA"].get_bfactor() if "CA" in r else 0.0)
        if ch.id == RECEPTOR_CHAIN:
            data["receptor_seq"] = "".join(seq)
            data["receptor_plddt"] = plddt
        elif ch.id == LIGAND_CHAIN:
            data["ligand_seq"] = "".join(seq)
            data["ligand_plddt"] = plddt
    return data


def compute_contact_metrics(structure, pae_matrix, idx_map):
    A = structure[0][RECEPTOR_CHAIN]
    B = structure[0][LIGAND_CHAIN]
    atomsA = [a for a in A.get_atoms() if a.element != "H"]
    atomsB = [a for a in B.get_atoms() if a.element != "H"]
    if not atomsA or not atomsB:
        return set(), [], ""
    ns = NeighborSearch(atomsA + atomsB)
    contacts = set()
    for a in atomsA:
        for b in ns.search(a.coord, CONTACT_DISTANCE_CUTOFF, level="A"):
            if b in atomsB:
                contacts.add((a.get_parent(), b.get_parent()))
    vals, lines = [], []
    for rA, rB in contacts:
        try:
            i = idx_map[(RECEPTOR_CHAIN, rA.id)]
            j = idx_map[(LIGAND_CHAIN, rB.id)]
        except KeyError:
            continue
        p = min(float(pae_matrix[i, j]), float(pae_matrix[j, i]))
        vals.append(p)
        lines.append(f"/{RECEPTOR_CHAIN}:{rA.id[1]} /{LIGAND_CHAIN}:{rB.id[1]} {p:.6f}")
    return contacts, vals, "\n".join(lines)


# ── Circular interaction plot ────────────────────────────────────────

def plot_interacting_residues(receptor_seq, ligand_seq, raw_contact_data_string,
                              output_filepath, title, **kwargs):
    def parse_contacts(raw):
        contacts, pae_scores = [], []
        pat = re.compile(r"^\s*/([A-Z]):(\d+)\s+/([A-Z]):(\d+)\s+(\d+\.\d+)$")
        for line in raw.strip().split("\n"):
            m = pat.match(line.strip())
            if not m:
                continue
            c1, r1, c2, r2, p = m.groups()
            r1, r2, p = int(r1), int(r2), float(p)
            if c1 == RECEPTOR_CHAIN and c2 == LIGAND_CHAIN:
                contacts.append((r1, r2))
                pae_scores.append(p)
            elif c1 == LIGAND_CHAIN and c2 == RECEPTOR_CHAIN:
                contacts.append((r2, r1))
                pae_scores.append(p)
        return contacts, pae_scores

    pairwise, pae_scores = parse_contacts(raw_contact_data_string)
    receptor_plddt = kwargs.get("receptor_plddt_scores", [])
    ligand_plddt = kwargs.get("ligand_plddt_scores", [])
    iptm_score = kwargs.get("iptm_score")
    avg_pae = kwargs.get("avg_pae")
    receptor_len = len(receptor_seq)
    ligand_len = len(ligand_seq)

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_aspect("equal")
    ax.set_axis_off()
    total_len = receptor_len + ligand_len
    gap_frac = 0.15
    rec_frac = receptor_len / total_len
    lig_frac = ligand_len / total_len
    rec_span = (1 - gap_frac) * rec_frac * 2 * np.pi
    lig_span = (1 - gap_frac) * lig_frac * 2 * np.pi
    gap = (gap_frac / 2) * 2 * np.pi
    rec_angles = np.linspace(gap, gap + rec_span, receptor_len)
    lig_angles = np.linspace(gap + rec_span + gap, 2 * np.pi, ligand_len)

    cmap_plddt = mcolors.ListedColormap(["#FF7D42", "#FFD100", "#65CBF3", "#0053D6"])
    norm_plddt = mcolors.BoundaryNorm([0, 50, 70, 90, 100], cmap_plddt.N)

    track_r, track_w = 1.0, 0.08
    for i, ang in enumerate(rec_angles):
        ax.add_patch(Wedge(
            (0, 0), track_r,
            np.rad2deg(ang - rec_span / receptor_len / 1.8),
            np.rad2deg(ang + rec_span / receptor_len / 1.8),
            width=track_w,
            facecolor=cmap_plddt(norm_plddt(receptor_plddt[i])),
            edgecolor="none",
        ))
    for i, ang in enumerate(lig_angles):
        ax.add_patch(Wedge(
            (0, 0), track_r,
            np.rad2deg(ang - lig_span / ligand_len / 1.8),
            np.rad2deg(ang + lig_span / ligand_len / 1.8),
            width=track_w,
            facecolor=cmap_plddt(norm_plddt(ligand_plddt[i])),
            edgecolor="none",
        ))

    pae_colors = ["#440154", "#3b528b", "#21918c", "#fde725", "#fd8d3c"]
    cmap_lines = mcolors.ListedColormap(pae_colors)
    # Fixed bins matching the legend: <2, 2-4, 4-6, 6-8, >8
    bounds = [0.0, 2.0, 4.0, 6.0, 8.0, 100.0]
    norm_lines = mcolors.BoundaryNorm(bounds, cmap_lines.N)
    for (r_res, l_res), p in zip(pairwise, pae_scores):
        a1 = rec_angles[r_res - 1]
        a2 = lig_angles[l_res - 1]
        x1, y1 = np.cos(a1) * (track_r - track_w), np.sin(a1) * (track_r - track_w)
        x2, y2 = np.cos(a2) * (track_r - track_w), np.sin(a2) * (track_r - track_w)
        path = MplPath(
            [(x1, y1), (0, 0), (x2, y2)],
            [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3],
        )
        ax.add_patch(PathPatch(
            path, facecolor="none",
            edgecolor=cmap_lines(norm_lines(p)),
            lw=1.8 if p < 2.0 else 0.5, alpha=0.9,
        ))

    total = len(pae_scores)
    under = sum(1 for p in pae_scores if p < SUMMARY_PAE_THRESHOLD)
    fig = plt.gcf()
    fig.text(
        0.5, 0.98,
        f"ipTM: {iptm_score if iptm_score is not None else float('nan'):.3f} | "
        f"Avg. PAE: {avg_pae:.2f} \u00c5 | Contacts < {SUMMARY_PAE_THRESHOLD}\u00c5: {under}/{total}",
        ha="center", va="top", fontsize=12, color="#333",
    )
    plddt_leg = [
        Rectangle((0, 0), 1, 1, fc=c, label=l) for c, l in
        zip(["#FF7D42", "#FFD100", "#65CBF3", "#0053D6"], ["<50", "50-70", "70-90", ">90"])
    ]
    fig.legend(handles=plddt_leg, title="pLDDT", loc="upper left",
               bbox_to_anchor=(0.02, 0.98), fontsize=9, title_fontsize=10)
    pae_leg = [
        Rectangle((0, 0), 1, 1, fc=c, label=l) for c, l in
        zip(pae_colors, ["<2\u00c5", "2-4\u00c5", "4-6\u00c5", "6-8\u00c5", ">8\u00c5"])
    ]
    fig.legend(handles=pae_leg, title="PAE", loc="upper right",
               bbox_to_anchor=(0.98, 0.98), fontsize=9, title_fontsize=10)
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    fig.text(0.5, 0.02, title, ha="center", va="bottom", fontsize=16)
    plt.savefig(output_filepath, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ── Analysis pipeline ────────────────────────────────────────────────

def compute_core_for_model0(src: ZipJobSource, job_name: str):
    cif_member = src.find(job_name, "*_model_0.cif")
    data_member = src.find(job_name, "*_full_data_0.json")
    if not (cif_member and data_member):
        raise FileNotFoundError(f"{job_name}: missing model_0 files")

    pae_matrix, full_json = load_pae(src, data_member)

    iptm = np.nan
    conf_member = src.find(job_name, "*_summary_confidences_0.json")
    if conf_member:
        try:
            iptm = float(json.loads(src.read_text(conf_member)).get("iptm", np.nan))
        except Exception:
            pass
    if np.isnan(iptm):
        iptm = extract_iptm(full_json)

    cif_text = src.read_text(cif_member)
    structure = parser_bio.get_structure(job_name, io.StringIO(cif_text))
    idx_map = residue_index_map(structure[0])
    contacts, pae_vals, raw = compute_contact_metrics(structure, pae_matrix, idx_map)
    avg_pae = float(np.mean(pae_vals)) if pae_vals else 0.0
    under = sum(1 for p in pae_vals if p < SUMMARY_PAE_THRESHOLD)
    total = len(pae_vals)

    return {
        "Job": job_name,
        "ipTM": float(iptm) if not np.isnan(iptm) else np.nan,
        "AvgPAE": avg_pae,
        "Contacts<5.0A PAE": f"{under} / {total}",
        "_structure": structure,
        "_raw": raw,
    }


def analyze_job_full(job_name: str, core: dict):
    structure = core["_structure"]
    seqs = get_sequences_and_plddt(structure)
    PLOTS_OUTPUT_DIR.mkdir(exist_ok=True)
    out_path = PLOTS_OUTPUT_DIR / f"{job_name}_interaction_plot.png"
    plot_interacting_residues(
        receptor_seq=seqs["receptor_seq"],
        ligand_seq=seqs["ligand_seq"],
        raw_contact_data_string=core["_raw"],
        output_filepath=out_path,
        title=f"Circular Interaction Summary for {job_name}",
        receptor_plddt_scores=seqs["receptor_plddt"],
        ligand_plddt_scores=seqs["ligand_plddt"],
        iptm_score=core["ipTM"],
        avg_pae=core["AvgPAE"],
    )
    print(f"  \u2713 {out_path}")
    return {
        "Job": core["Job"],
        "ipTM": core["ipTM"],
        "AvgPAE": core["AvgPAE"],
        "Contacts<5.0A PAE": core["Contacts<5.0A PAE"],
    }


# ── Main ─────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="AlphaScan — AlphaFold3 result analysis from zip archives")
    ap.add_argument("zip_files", nargs="*", type=Path,
                    help="AF3 zip file(s). If omitted, uses all .zip in current directory.")
    ap.add_argument("--iptm-cutoff", type=float, default=0.0,
                    help="Keep jobs with ipTM >= this value (default: 0.0)")
    ap.add_argument("--list-iptm", action="store_true",
                    help="Quick mode: show top ipTM scores and exit.")
    ap.add_argument("--top", type=int, default=20,
                    help="Number of rows in quick list (default: 20)")
    args = ap.parse_args()
    cutoff = args.iptm_cutoff

    # Resolve zip files
    zip_paths = args.zip_files
    if not zip_paths:
        zip_paths = sorted(Path(".").glob("*.zip"))
        if not zip_paths:
            print("No .zip files found in the current directory.")
            sys.exit(1)
        print(f"[info] Found {len(zip_paths)} zip file(s) in current directory")

    sources: list[ZipJobSource] = []
    for zp in zip_paths:
        if not zp.exists():
            print(f"[warn] zip not found, skipping: {zp}")
            continue
        src = ZipJobSource(zp)
        jn = [j for j in src.job_names() if not j.startswith(".")]
        if jn:
            sources.append(src)
            print(f"  {zp.name}: {len(jn)} jobs")
        else:
            print(f"  {zp.name}: no job folders, skipping")
            src.close()

    if not sources:
        print("No job folders found in any zip.")
        return

    all_jobs = []
    for src in sources:
        for jn in src.job_names():
            if not jn.startswith("."):
                all_jobs.append((src, jn))

    # Quick list mode
    if args.list_iptm:
        rows = []
        for src, jn in all_jobs:
            try:
                s = read_model0_iptm_fast(src, jn)
                if not np.isnan(s):
                    rows.append({"Job": jn, "ipTM": float(s)})
            except Exception as e:
                print(f"[SKIP] {jn}: {e}")
        if not rows:
            print("No ipTM values found.")
        else:
            df = pd.DataFrame(rows).sort_values("ipTM", ascending=False).head(args.top)
            with pd.option_context("display.float_format", "{:.4f}".format):
                print(df.to_string(index=False))
        for src in sources:
            src.close()
        return

    # Full scan
    kept = []
    for src, jn in all_jobs:
        try:
            iptm_fast = read_model0_iptm_fast(src, jn)
            if np.isnan(iptm_fast):
                if cutoff > 0:
                    continue
            elif iptm_fast < cutoff:
                continue
            kept.append((src, jn, iptm_fast))
        except Exception as e:
            print(f"[SKIP] {jn}: {e}")

    if not kept:
        print(f"No jobs meet ipTM cutoff ({cutoff}).")
        for src in sources:
            src.close()
        return

    rows = []
    cores = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        for src, jn, _ in kept:
            try:
                core = compute_core_for_model0(src, jn)
                rows.append({k: core[k] for k in ("Job", "ipTM", "AvgPAE", "Contacts<5.0A PAE")})
                cores[jn] = core
            except Exception as e:
                print(f"[SKIP] {jn}: {e}")

    if not rows:
        print("No jobs to summarise.")
        for src in sources:
            src.close()
        return

    df = pd.DataFrame(rows)
    df = df.sort_values(by=["ipTM", "AvgPAE"], ascending=[False, True], na_position="last")
    with pd.option_context("display.float_format", "{:.4f}".format):
        print(df[["Job", "ipTM", "AvgPAE", "Contacts<5.0A PAE"]].to_string(index=False))

    max_n = len(df)
    ans = input(f"\nHow many to fully analyse (plots + CSV)? [0-{max_n}] (Enter=all, 0=none): ").strip()
    n = max_n if ans == "" else (int(ans) if ans.isdigit() else 0)
    n = max(0, min(n, max_n))
    if n == 0:
        print("\nNo jobs selected for full analysis.")
        for src in sources:
            src.close()
        return

    results = []
    selected = set(df["Job"].head(n))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        for dname in selected:
            try:
                results.append(analyze_job_full(dname, cores[dname]))
            except Exception as e:
                print(f"  [FAIL] {dname}: {e}")

    if results:
        pd.DataFrame(results)[["Job", "ipTM", "AvgPAE", "Contacts<5.0A PAE"]] \
            .sort_values("ipTM", ascending=False, na_position="last") \
            .to_csv(CORE_SUMMARY_FILENAME, index=False, float_format="%.4f")
        print(f"\n\u2713 Saved: {CORE_SUMMARY_FILENAME}")

    for src in sources:
        src.close()


if __name__ == "__main__":
    main()
