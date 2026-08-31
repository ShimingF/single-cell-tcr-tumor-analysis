# Data

## Source Dataset

This project uses the publicly available **NCBI Gene Expression Omnibus dataset GSE139555**:

> *Peripheral clonal expansion of T lymphocytes associates with tumour infiltration and response to cancer immunotherapy*

NCBI GEO accession: **GSE139555**  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139555

The GEO series contains paired **single-cell RNA-seq (scRNA-seq)** and **single-cell TCR-seq (scTCR-seq)** data from pretreatment samples of 14 cancer patients. Samples include tumor, normal-adjacent tissue (NAT), and peripheral blood where available.

Processed supplementary data are available through GEO. According to the GEO record, raw data were deposited in the European Genome-phenome Archive (EGA):

- scRNA-seq: **EGAS00001003993**
- scTCR-seq: **EGAS00001003994**

## Why the Data Are Not Included Here

Raw and processed source datasets are intentionally **not redistributed in this GitHub repository**.

Reasons include:

- large file sizes;
- preservation of provenance to the original public repositories;
- avoidance of unnecessary redistribution of patient-derived molecular data;
- keeping the repository focused on analysis code and reproducible workflow structure.

## Expected Local Structure

A local analysis environment may use a structure similar to:

```text
data/
├── raw/
│   └── TCR/
└── processed/
    ├── final_paired_tcell_metadata.rds
    └── tcrex_results_score_gt_0.5.csv
```

These files and directories are excluded from Git version control by `.gitignore`.

## Reproduction

1. Download the appropriate GSE139555 supplementary files from NCBI GEO.
2. Obtain any controlled raw data from EGA if required and if access is authorized.
3. Place the source files in the local `data/` directory expected by the analysis scripts.
4. Run the preprocessing and analysis scripts in numerical order.

The derived file `final_paired_tcell_metadata.rds` is generated locally by the analysis pipeline and should not be committed to GitHub.

## Attribution

All original GSE139555 data remain the work of the original dataset contributors. This repository contains an independent downstream computational reanalysis and does not claim ownership of the source data.
