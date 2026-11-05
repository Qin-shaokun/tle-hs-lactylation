# Data

No expression data are stored in this repository. Everything is public; the table below gives the source, the access route and what the pipeline expects to find on disk.

## Accessions

| Accession | Study | Samples used here | Route |
|---|---|---|---|
| GSE256068 | François et al., *Nat Commun* 2024;15:2180 | 64 TLE-HS + 14 non-epileptic controls, cortex | `GEOquery::getGEO("GSE256068")`; processed matrix from `suppl/` |
| GSE190451 | Chen et al., *Nat Neurosci* 2023;26:542–554 | 3 TLE + 3 controls, cortex | `GEOquery::getGEO("GSE190451")` |
| GSE71058 | Griffin et al., *Epilepsia* 2016;57:376–385 | Laser-capture dentate gyrus granule cells, HS vs no-HS | `GEOquery::getGEO("GSE71058")` |
| GSE190452 | Chen et al., *Nat Neurosci* 2023;26:542–554 | Single-nucleus RNA-seq, temporal cortex, 4 TLE + 4 controls | Raw archive: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190452/suppl/GSE190452_RAW.tar` (~698 MB) |

GSE190451 and GSE190452 are two sub-series of the same study (bulk and single-nucleus respectively).

## GSE190452 sample map

The 4 + 4 design was submitted in two batches, which matters for the QC and clustering steps:

| GSM | Label | Condition | Batch |
|---|---|---|---|
| GSM5724021 | Normal_1 | Control | 1 |
| GSM5724022 | Normal_2 | Control | 1 |
| GSM5724023 | Normal_3 | Control | 1 |
| GSM5724024 | TLE_1 | TLE | 1 |
| GSM5724025 | TLE_2 | TLE | 1 |
| GSM5724026 | TLE_3 | TLE | 1 |
| GSM6893724 | Normal_4 | Control | 2 |
| GSM6893725 | TLE_4 | TLE | 2 |

Condition is balanced across batches, so batch is confounded with neither group nor with the four-versus-four split as a whole — but it is not ignorable, and the UMAP is coloured by batch so that the reader can judge.

## Expected layout after download

```
data/
├── GSE256068/            # created by R/01_download_and_metadata.R
├── GSE190451/
├── GSE71058/
├── GSE190452_raw/        # unpacked 10x triplets, one flat directory
│   ├── GSM5724021_Normal1_barcodes.tsv.gz
│   ├── GSM5724021_Normal1_features.tsv.gz
│   ├── GSM5724021_Normal1_matrix.mtx.gz
│   └── ... (8 samples × 3 files)
└── lactylation_geneset_v1.csv
```

## Notes on the gene set

`lactylation_geneset_v1.csv` lists 49 genes in six categories (writers, erasers, readers, lactate metabolism, lactate transport, substrate genes). Each row carries a `review_status` field; entries marked `pending` were retained on the basis of the literature but not individually re-verified against a primary source. The heat map reports the 47 genes detectable in GSE256068.

One caveat worth recording: `RFWD2` appears as the alias **COP1** in the GEO feature annotation, and the feature table contains 10 duplicated gene symbols (two Ensembl IDs per symbol), which `R/09_sc_qc.R` resolves with `make.unique()` before constructing the Seurat object.
