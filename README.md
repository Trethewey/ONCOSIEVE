<div align="center">
  <h1>O N C O S I E V E</h1>
  <p>Pan-cancer variant curation and rescue tool</p>
  <p><em>7 databases · 33.5 million variants · 30,547 curated whitelist entries</em></p>
  <img src="https://github.com/user-attachments/assets/56304b21-8193-4e84-869a-347aadf7ab76" width="450"/>
  <p><strong>Author:</strong> Dr Christopher Trethewey<br>
  <strong>Email:</strong> christopher.trethewey@nhs.net</p>
</div>

---

## Overview

ONCOSIEVE builds a pan-cancer somatic variant whitelist from multiple curated
databases and applies it to rescue clinically relevant low-VAF variants from
Mutect2 post-filter calls. Designed for use in research NGS pipelines.

> [!WARNING]
> **Research use only.** ONCOSIEVE is an experimental research tool and has not
> been validated for clinical diagnostic use. It must not be used to inform
> patient management or clinical decision-making without independent validation
> in an accredited diagnostic setting.

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

Each variant in the whitelist is assigned a confidence tier based on the
evidence supporting it. Tiering determines the minimum VAF required to rescue
a Mutect2-filtered variant: higher confidence variants can be rescued at lower
VAFs.

A variant is assigned the highest tier for which it qualifies.

| Tier | Criteria | Min VAF for rescue |
|------|----------|--------------------|
| 1 | OncoKB Oncogenic or Likely Oncogenic, OR seen in ≥50 samples across ≥3 cancer types | 0.5% |
| 2 | OncoKB Predicted Oncogenic, OR present in ClinVar (pathogenic/somatic), OR a CancerHotspots recurrent mutation, OR seen in ≥25 samples across ≥2 cancer types | 0.5% |
| 3 | Count-based only, below Tier 2 count threshold | 1.0% |

**Why Tier 3 is currently empty:** The base filter for inclusion in the whitelist
requires n≥25 samples and ct≥2 cancer types. The Tier 2 count threshold is set
to the same values (n≥25, ct≥2). Any variant that passes the base filter
therefore automatically meets the Tier 2 count threshold, leaving nothing to
fall through to Tier 3. Tier 3 would populate if the base filter were relaxed
below the Tier 2 threshold in `settings.yaml`.

VAF floors per tier are configurable in `settings.yaml` under `vaf_rescue`.

---

## Transcript annotation

HGVSc strings in the whitelist are carried through from each source database
via VEP annotation. All coordinate-resolved variants are annotated against
MANE Select transcripts during the VEP remapping steps in data preparation
(see `tools/hotspots_vep_remap.py`). HGVSc values use ENST accessions
throughout.

The `transcript_id` column in the annotated TSV contains the ENST accession
extracted from `hgvsc`. The `protein_change` column contains the normalised
3-letter HGVS protein change derived from `hgvsp`.

**Note on `tools/mane_remap.py`:** This tool is not currently compatible with
the pipeline output. It expects RefSeq NM_ accessions but the pipeline produces
ENST-prefixed HGVSc values from VEP. It is retained for potential future use
but should not be run as a post-processing step.

---

## Configuration

Edit `settings.yaml` to change thresholds, VAF floors, and OncoKB token.
Edit `config.yaml` to change file paths or disable sources.

---

## OncoKB token

Academic licence. Token expires in 6 months. Update in `settings.yaml` under `oncokb.api_token`.
