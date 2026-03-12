<div align="center">
  <h1>O N C O S I E V E</h1>
  <p><em>Pan-cancer variant curation and rescue tool</em></p>
  <img src="https://github.com/user-attachments/assets/56304b21-8193-4e84-869a-347aadf7ab76" width="450"/>
  <p><strong>Author:</strong> Dr Christopher Trethewey<br>
  <strong>Email:</strong> christopher.trethewey@nhs.net</p>
</div>

---

## Overview

ONCOSIEVE builds a pan-cancer somatic variant whitelist from multiple curated
databases and applies it to rescue clinically relevant low-VAF variants from
Mutect2 post-filter calls. Designed for use in diagnostic NGS pipelines.

---

## Data sources

| Source          | Version   | Type | Approx. variants (source total) | Notes                                                                  |
|-----------------|-----------|------|---------------------------------|------------------------------------------------------------------------|
| COSMIC          | v103      | File | ~25,000,000                     | GRCh38 TSV + VCF; >6,800 cancer types                                 |
| AACR GENIE      | v19.0     | File | ~3,750,000                      | MAF; 271,837 samples / 227,696 patients; 19 cancer centres            |
| ClinVar         | 2025      | File | ~1,000,000                      | GRCh38 VCF; pathogenic/likely pathogenic and somatic-flagged filtered |
| OncoKB          | Current   | API  | ~7,700                          | ~850 genes; 130 cancer types; academic token required                  |
| TP53 database   | R21       | File | ~29,900                         | GRCh38 CSV; functional annotations for >9,000 mutant proteins         |
| cBioPortal      | Current   | API  | ~3,750,000*                     | Live REST API; 532 hg38-compatible studies queried                    |
| CancerHotspots  | v2        | API  | ~3,181                          | Live REST API; 24,592 tumour samples; q-value filtered                |
|                 |           |      |                                 |                                                                        |
| **Raw total**   |           |      | **~33,540,781**                 | Pre-deduplication; significant inter-database overlap expected         |

`*` cBioPortal figure reflects mutations retrieved by this pipeline, not the portal's full dataset.

---

## Directory structure
```
oncosieve/
├── build_whitelist.py       # Main pipeline entry point
├── mutect2_rescue.py        # Mutect2 post-filter rescue
├── pre_check.py             # Pre-run dependency and data audit
├── test_pipeline.py         # Test harness
├── run_pipeline.sh          # Orchestrator; run this to start the pipeline
├── config.yaml              # File paths and source enable/disable flags
├── settings.yaml            # Thresholds and parameters (edit this)
├── requirements.txt
├── panels/
│   └── lymphoma.bed         # Add panel BED files here
├── parsers/
│   ├── parse_cosmic.py
│   ├── parse_genie.py
│   ├── parse_clinvar.py
│   ├── parse_oncokb.py
│   ├── parse_tp53.py
│   ├── parse_cbioportal.py
│   └── parse_hotspots.py
├── tools/                   # Standalone utilities; run manually
│   ├── annotate_panels.py
│   ├── clinvar_vep_annotate.py
│   ├── hotspots_vep_remap.py
│   ├── mane_audit.py
│   └── mane_remap.py
└── logs/                    # Auto-created on first run
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
python tools/annotate_panels.py \
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

| Tier | Criteria                                                          | Min VAF |
|------|-------------------------------------------------------------------|---------|
| 1    | OncoKB Oncogenic/Likely Oncogenic, or n≥50, ct≥3                 | 0.5%    |
| 2    | OncoKB Predicted Oncogenic, ClinVar, Hotspots                    | 0.5%    |
| 3    | Count-based only (n≥25, ct≥2)                                    | 1.0%    |

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
python tools/mane_remap.py \
    --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \
    --output    output/pan_cancer_whitelist_GRCh38.mane.tsv.gz

# Test on first 500 variants
python tools/mane_remap.py \
    --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \
    --output    output/pan_cancer_whitelist_GRCh38.mane.tsv.gz \
    --max-variants 500
```

The output adds two columns: `hgvsc_mane` (MANE Select HGVSc where available,
original otherwise) and `mane_remapped` (True/False). Variants where remapping
fails via the API retain their original HGVSc value.

Requires: `pip install requests pandas`

---

## Configuration

Edit `settings.yaml` to change thresholds, VAF floors, and OncoKB token.
Edit `config.yaml` to change file paths or disable sources.

---

## OncoKB token

Academic licence. Token expires in 6 months. Update in `settings.yaml` under `oncokb.api_token`.
