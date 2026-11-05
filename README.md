# Glycolytic–lactylation programme in human TLE-HS cortex

Analysis code for the study *"A glycolytic–lactylation transcriptional programme in human temporal lobe epilepsy with hippocampal sclerosis"*.

The repository contains every script needed to reproduce the analyses, figures and supplementary figures from publicly available data. Raw and processed expression data are **not** included; the download step is scripted (see [Data](#data)).

---

## What the study does

Bulk RNA-sequencing of 78 cortical specimens (64 TLE-HS, 14 non-epileptic controls; GSE256068) is used to ask three questions:

1. Does the epileptogenic cortex carry a coordinated glycolytic–lactylation transcriptional signature?
2. Which cell compartment does it reside in?
3. Does it coexist with suppression of programmes required for cognition (synaptic vesicle machinery, ribosomal biogenesis)?

The pipeline covers covariate-adjusted differential expression, WGCNA, consensus machine-learning feature selection, Firth-penalised nomogram with bootstrap calibration, ssGSEA cell-type deconvolution, ORA/GSEA enrichment, and single-nucleus localisation in an independent cohort (GSE190452).

---

## Data

All datasets are public. `R/01_download_and_metadata.R` retrieves them via `GEOquery`; large raw archives are downloaded separately (see `data/README.md`).

| Dataset | Role | Content |
|---|---|---|
| **GSE256068** | Primary discovery | Human cortex, 64 TLE-HS + 14 controls (François et al., *Nat Commun* 2024) |
| **GSE190451** | Independent replication | Human cortex, 3 TLE + 3 controls (Chen et al., *Nat Neurosci* 2023) |
| **GSE71058** | Severity replication | Laser-capture dentate gyrus granule cells, HS vs no-HS (Griffin et al., *Epilepsia* 2016) |
| **GSE190452** | Cell-type localisation | Single-nucleus RNA-seq, temporal cortex, 4 TLE + 4 controls (Chen et al., *Nat Neurosci* 2023) |

---

## Requirements

- R ≥ 4.6 (developed on 4.6.2; exact package versions in `docs/session_info.txt`)
- CRAN: `glmnet`, `e1071`, `randomForest`, `pROC`, `rms`, `logistf`, `dplyr`, `ggplot2`, `plyr`, `patchwork`
- Bioconductor: `GEOquery`, `limma`, `WGCNA`, `GSVA`, `GSEABase`, `clusterProfiler`, `org.Hs.eg.db`, `Seurat`

```r
install.packages(c("glmnet","e1071","randomForest","pROC","rms","logistf",
                   "dplyr","ggplot2","plyr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery","limma","WGCNA","GSVA","GSEABase",
                       "clusterProfiler","org.Hs.eg.db","Seurat"))
```

---

## Running

```bash
git clone <repository-url>
cd <repository>
Rscript run_all.R          # full pipeline, fails fast
```

`run_all.R` checks that all packages are present before starting, then sources each script in order. Working directory must be the repository root — paths in `config.R` are relative to it. All thresholds and dataset parameters live in `config.R`; the scripts themselves contain no hard-coded paths.

Individual steps can be run in isolation (`Rscript R/02_deg_main.R`) as long as the preceding step has written its intermediate object to `results/`.

### Pipeline order

| Script | Purpose | Key output |
|---|---|---|
| `01_download_and_metadata.R` | Fetch GEO series, parse sample metadata | `results/01_*` |
| `02_deg_main.R` | limma, `~ age + sex + condition` | `02_deg_all.csv`, volcano, heatmap |
| `02b_sensitivity.R` | Adult subset and region-adjusted models | `S1_sensitivity.png` |
| `03_wgcna.R` | Signed network, power 8, 5,000 most variable genes | module–trait heatmap |
| `04_ml_hub_selection.R` | LASSO, random forest, SVM-RFE | selection tables, ROC |
| `04b_ml_venn_rf.R` | Three-way overlap and RF importance ranking | Venn, importance plot |
| `05_validation_roc.R` | ROC in the replication cohort | `05_roc_valid.png`, `S5_roc_severity.png` |
| `05b_nomogram_firth.R` | Firth logistic model, nomogram | `05_nomogram.png` |
| `05c_calibration.R` | Bootstrap calibration (Harrell) | `05_calibration.png` |
| `06_immune_infiltration.R` | ssGSEA cell-type scores, gene–cell correlation | boxplots, correlation heatmap |
| `08_enrichment.R` | GO/KEGG ORA and GSEA | enrichment and GSEA panels |
| `08b_kegg.R` | KEGG GSEA with retry logic | `08_gsea_KEGG.png` |
| `09_sc_qc.R` | snRNA-seq QC (GSE190452) | per-sample QC tables, violin plots |
| `09b_cluster.R` | PCA/UMAP, clustering, marker-based annotation | annotated Seurat object |
| `09b3_composition_hub.R` | Per-sample composition, hub-gene detection rate | composition panels |
| `09c_localize.R` | Cell-type localisation of hub and lactylation genes | localisation tables |
| `09c_keygenes_dotplot.R` | Key-gene dot plot | dot plot |
| `09c2_markers_fast.R` | Marker check on variable features (annotation QC) | marker table |
| `09d_assemble_figs.R`, `09e_print_panels.R` | Panel assembly, print-resolution redraw | Figure 7, S6, S7 |
| `11_redraw_panels.R` | Redraw panels at final print size | print panels |
| `10_assemble_main.R` | Assemble Figures 1–6 at 180 mm / 300 dpi | `figures/assembled/` |

---

## Analytical decisions worth knowing

These are deliberate, and each is stated in the paper's Methods:

- **Design matrix** is `~ age + sex + condition`; two prespecified sensitivity analyses add sampling region and restrict to adults (≥18 years). The control group contains 10 donors under 18, so age adjustment is not cosmetic.
- **Hub genes** are the intersection of LASSO and random-forest selections. The three-algorithm intersection is empty and is shown as such (Figure 3C); SVM-RFE is reported but does not define the markers.
- **Cell-type signatures exclude the hub genes**, to avoid correlating a gene with a signature that contains it. Omitting this step inflates the B-cell correlation from ρ = 0.33 to ρ = 0.81.
- **Microglial nuclei are scarce in control specimens (n = 7)**, so no within-microglia group comparison is performed; the single-nucleus data localise the signal but do not establish per-cell-type differential expression.
- **Cell-composition estimates are descriptive only.** Some cell types are absent from individual specimens, so group means are reported alongside the per-specimen range, without testing.
- **The glycolytic and lactylation signal is transcriptional.** No lactylome, no pan-Kla immunoblot and no cognitive scores are available; the paper explicitly frames the result as a testable transcriptional framework rather than a biochemical demonstration.

---

## Repository layout

```
.
├── config.R              # all paths, dataset parameters, thresholds
├── run_all.R             # pipeline driver
├── R/                    # analysis scripts, numbered by execution order
├── data/                 # gene set table + download instructions (no raw data)
├── docs/                 # session info, package versions
├── LICENSE
└── README.md
```

Outputs are written to `results/` (tables, `.rds` objects) and `figures/` (`figures/assembled/` holds the submission-resolution versions).

---

## Data and code availability

Expression data are available from GEO under the accessions listed above. This repository contains all custom code; no custom algorithms were developed. Figures were assembled with base R `grid` and `png`.

## License

MIT — see `LICENSE`.

## Citation

If you use this code, please cite the paper (see `CITATION.cff`) and the original data papers listed in `data/README.md`.
