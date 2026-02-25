# 🔬 AlphaScan

**AlphaFold3 Result Analysis Pipeline** — reads directly from AF3 zip archives, no extraction needed.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)

---

AlphaScan is a command-line tool that takes the zip files you download from the [AlphaFold3 Server](https://alphafoldserver.com/) and automatically ranks, scores, and visualises your protein interaction predictions. It saves you from manually opening dozens of result folders — just point it at your zips and get a ranked summary with publication-ready circular interaction plots. 🎯

## ✨ Features

- 🚀 **Quick ipTM scan** — fast ranking of all jobs by ipTM score (interface predicted TM-score)
- 📊 **Full analysis** — ipTM, average PAE, contact metrics, and circular interaction plots
- 📦 **Reads zip archives directly** — no need to extract AF3 results first
- 🗂️ **Handles flat and nested zips** — auto-detects job folder structure
- 🎛️ **Interactive top-N selection** — scan all jobs, then choose how many to fully analyse

## 🧬 What Does It Measure?

AlphaScan extracts three key metrics from AlphaFold3 predictions:

| Metric | What It Means | Why It Matters |
|---|---|---|
| **ipTM** | Interface predicted TM-score (0–1) | Higher = more confident that the predicted protein–protein interface is correct. Scores above **0.5** are generally promising |
| **PAE** | Predicted Aligned Error (Å) | Lower = AlphaFold is more confident about the relative positions of residues. Contact lines with PAE < 5Å are high-confidence |
| **pLDDT** | Per-residue confidence (0–100) | Higher = more confident local structure. Blue (>90) is very high, orange (<50) is low |

## 🗺️ Pipeline Workflow

```
📂 Point at AF3 zip files
     │
     ▼
🚀  Quick ipTM Scan  (reads scores from all jobs — fast!)
     │
     ▼
🔍  Filter  (optional: keep only jobs above your ipTM cutoff)
     │
     ▼
📋  Ranking Table  (sorted by ipTM with PAE and contact counts)
     │
     ▼
🎛️  Select Top-N  (choose how many top jobs to fully analyse)
     │
     ▼
🎨  Generate Plots  (circular interaction plots for selected jobs)
     │
     ▼
💾  Save Results  (summary CSV + plot PNGs)
```

### Step-by-step

1. **Quick ipTM scan** — AlphaScan opens each zip file and reads the `model_0` ipTM score from the JSON confidence data. This is fast because it skips the heavy structure files
2. **Filtering** — Optionally discard jobs below a user-defined ipTM cutoff (e.g. `--iptm-cutoff 0.5` to keep only promising hits)
3. **Ranking table** — All qualifying jobs are sorted by ipTM (highest first) and displayed in a table with ipTM, average PAE, and the number of inter-chain contacts with PAE below 5.0 Å
4. **Interactive prompt** — You choose how many of the top-ranked jobs to fully analyse (e.g. "analyse top 10")
5. **Full analysis** — For each selected job, AlphaScan parses the CIF structure file, identifies inter-chain atomic contacts within 4.5 Å, and generates a circular interaction plot (circos-style)

## 📋 Requirements

- Python 3.8 or later
- The following Python packages:

```bash
pip install numpy pandas matplotlib biopython
```

## 🚀 Usage

```bash
# Analyse all .zip files in current directory
python alphascan.py

# Analyse specific zip files
python alphascan.py results1.zip results2.zip

# Quick ipTM ranking only (top 30)
python alphascan.py --list-iptm --top 30

# Filter by ipTM cutoff (keep only promising jobs)
python alphascan.py --iptm-cutoff 0.5
```

### ⚙️ Arguments

| Argument | Description | Default |
|---|---|---|
| `zip_files` | AF3 zip file(s). If omitted, uses all `.zip` in current directory | All `.zip` |
| `--iptm-cutoff` | Keep jobs with ipTM ≥ this value | `0.0` |
| `--list-iptm` | Quick mode: show top ipTM scores and exit | Off |
| `--top` | Number of rows in quick list | `20` |

## 📤 Output

AlphaScan creates two outputs in the same directory as your zip files:

- 🎨 **`final_interaction_plots/`** — circular interaction plot PNGs for each analysed job
- 📄 **`afm_core_scores.csv`** — summary CSV with ipTM, average PAE, and contact counts

### 🎨 Circular Interaction Plot

Each plot is a circos-style visualisation of how two protein chains interact:

- **Outer ring** — per-residue pLDDT confidence, coloured by quality:
  - 🟠 Orange: < 50 (low confidence)
  - 🟡 Yellow: 50–70 (moderate)
  - 🔵 Cyan: 70–90 (good)
  - 🔷 Blue: > 90 (very high confidence)
- **Inner lines** — inter-chain contacts coloured by PAE (predicted error):
  - 🟣 Purple: < 2 Å (excellent)
  - 🔵 Blue: 2–4 Å (good)
  - 🟢 Teal: 4–6 Å (moderate)
  - 🟡 Yellow: 6–8 Å (uncertain)
  - 🟠 Orange: > 8 Å (low confidence)
- **Header** — ipTM score, average PAE, and number of contacts below threshold

## 👥 Authors

- **Nick Eilmann** — 📧 [nme122@ic.ac.uk](mailto:nme122@ic.ac.uk) · 🐙 [@ne1899](https://github.com/ne1899)
- **Tolga Bozkurt** — 📧 [o.bozkurt@ic.ac.uk](mailto:o.bozkurt@ic.ac.uk)

## 📝 Citation

If you use AlphaScan in your research, please cite:

> Eilmann N., Bozkurt T. (2026). AlphaScan: AlphaFold3 Result Analysis Pipeline. https://github.com/ne1899/Alphascan

## 📄 License

This project is licensed under the MIT License &mdash; see the [LICENSE](LICENSE) file for details.
