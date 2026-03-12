# <div align="center">
  <h1>O N C O S I E V E</h1>
  <p><em>Pan-cancer variant curation and rescue tool</em></p>
</div>

## Pan-cancer variant curation and rescue tool

![9c87c033-dd0d-42fb-a4d5-616da12bbed8](https://github.com/user-attachments/assets/56304b21-8193-4e84-869a-347aadf7ab76)

**Author:** Dr Christopher Trethewey  
**Email:** christopher.trethewey@nhs.net

---

## Overview

ONCOSIEVE builds a pan-cancer somatic variant whitelist from multiple curated
databases and applies it to rescue clinically relevant low-VAF variants from
Mutect2 post-filter calls. Designed for use in NHS diagnostic NGS pipelines.

---

## Data sources

| Source         | Type        | Notes                                      |
|----------------|-------------|--------------------------------------------|
| COSMIC v103    | File        | GRCh38 TSV + VCF                           |
| AACR GENIE v19 | File        | MAF, 271,837 samples                       |
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


---

## Transcript annotation and MANE Select

HGVSc strings in the whitelist are carried through as-is from each source
database. COSMIC, GENIE, TP53, and cBioPortal each use their own reference
transcripts, which may differ between sources and may not correspond to the
MANE Select transcript for a given gene.

**Known limitation:** No transcript normalisation is applied during the main
pipeline. Two entries for the same variant from different sources may carry
different HGVSc strings if those sources used different transcripts.

A post-processing script is provided to remap HGVSc annotations to MANE Select
transcripts using the NCBI MANE Select table and the Ensembl VEP REST API:

```bash
python mane_remap.py \
    --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \
    --output    output/pan_cancer_whitelist_GRCh38.mane.tsv.gz

# Test on first 500 variants
python mane_remap.py \
    --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \
    --output    output/pan_cancer_whitelist_GRCh38.mane.tsv.gz \
    --max-variants 500
```

The output adds two columns: `hgvsc_mane` (MANE Select HGVSc where available,
original otherwise) and `mane_remapped` (True/False). Variants where remapping
fails via the API retain their original HGVSc value.

Requires: `pip install requests pandas`

## Configuration

Edit `settings.yaml` to change thresholds, VAF floors, and OncoKB token.  
Edit `config.yaml` to change file paths or disable sources.

---

## OncoKB token

Academic licence. Token expires in 6 months. Update in `settings.yaml` under `oncokb.api_token`.
