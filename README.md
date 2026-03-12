# ONCOSIEVE
## Pan-cancer variant curation and rescue tool

**Author:** Dr Christopher Trethewey  
**Email:** christopher.trethewey@nhs.net

---

## Overview

ONCOSIEVE builds a pan-cancer somatic variant whitelist from multiple curated
databases and applies it to rescue clinically relevant low-VAF variants from
Mutect2 post-filter calls. Designed for use in future diagnostic NGS pipelines.

---

## Data sources

| Source         | Type        | Notes                                      |
|----------------|-------------|--------------------------------------------|
| COSMIC v103    | File        | GRCh38 TSV + VCF                           |
| AACR GENIE     | File        | MAF, 271,837 samples                       |
| ClinVar        | File        | GRCh38 VCF                                 |
| OncoKB         | API         | Annotation endpoint; academic token        |
| TP53 database  | File        | GRCh38 CSV                                 |
| cBioPortal     | API         | Live REST API                              |
| CancerHotspots | API         | Live REST API                              |
| DoCM           | API         | hg19 coordinates; flagged accordingly      |

---

## Directory structure

```
oncosieve/
├── build_whitelist.py       # Main pipeline
├── mutect2_rescue.py        # Mutect2 post-filter rescue
├── annotate_panels.py       # Panel membership annotation
├── run_pipeline.sh          # Entry point
├── config.yaml              # File paths and source toggles
├── settings.yaml            # Thresholds and parameters (edit this)
├── requirements.txt
├── panels/
│   └── lymphoma.bed         # Add panel BED files here
└── parsers/
    ├── parse_cosmic.py
    ├── parse_genie.py
    ├── parse_clinvar.py
    ├── parse_oncokb.py
    ├── parse_tp53.py
    ├── parse_cbioportal.py
    ├── parse_hotspots.py
    └── parse_docm.py
```

---

## Usage

### 1. Build the whitelist

```bash
source ~/venv_ngs/bin/activate
cd /path/to/oncosieve
bash run_pipeline.sh /path/to/data/
```

### 2. Annotate with diagnostic panels

```bash
python annotate_panels.py \
    --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \
    --panels-dir panels/ \
    --output output/pan_cancer_whitelist_GRCh38.annotated.tsv.gz
```

Add panels by dropping BED files into `panels/`. The script picks them up automatically.

### 3. Rescue Mutect2 variants

```bash
python mutect2_rescue.py \
    --input sample.vcf.gz \
    --whitelist output/pan_cancer_whitelist_GRCh38.vcf.gz \
    --output sample.rescued.vcf.gz
```

---

## Tiering

| Tier | Criteria                                                              | Min VAF |
|------|-----------------------------------------------------------------------|---------|
| 1    | OncoKB Oncogenic/Likely Oncogenic, or n≥50, ct≥3                     | 0.5%    |
| 2    | OncoKB Predicted Oncogenic, ClinVar, Hotspots, DoCM                  | 0.5%    |
| 3    | Count-based only (n≥25, ct≥2)                                        | 1.0%    |

---

## Configuration

Edit `settings.yaml` to change thresholds, VAF floors, and OncoKB token.  
Edit `config.yaml` to change file paths or disable sources.

---

## OncoKB token

Academic licence. Token expires in 6 months. Update in `settings.yaml` under `oncokb.api_token`.
