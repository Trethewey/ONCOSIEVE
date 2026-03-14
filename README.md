<div align="center">
  <h1>ONCOSIEVE</h1>
  <p><strong>Pan-cancer somatic variant whitelist for Mutect2 post-filter rescue</strong></p>
  <p>7 databases · 33,391,744 source variants · 33,971 curated whitelist entries</p>
  <img src="https://github.com/user-attachments/assets/56304b21-8193-4e84-869a-347aadf7ab76" width="450"/>
  <p>
    <strong>Author:</strong> Dr Christopher Trethewey<br>
    <strong>Email:</strong> christopher.trethewey@nhs.net
  </p>
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

| Source          | Version      | Type | Approx. variants (source total) | Notes                                                                              |
|-----------------|--------------|------|---------------------------------|------------------------------------------------------------------------------------|
| COSMIC          | v103         | File | ~25,000,000                     | GRCh38 TSV + VCF; >6,800 cancer types                                             |
| AACR GENIE      | v19.0        | File | ~3,750,000                      | MAF; 271,837 samples / 227,696 patients; 19 cancer centres; GRCh37 lifted to GRCh38 |
| TCGA mc3        | v0.2.8       | File | ~3,600,963                      | PanCancer Atlas MAF; 33 cancer types; 10,295 tumours; GRCh37 lifted to GRCh38     |
| ClinVar         | 2025         | File | ~1,000,000                      | GRCh38 VCF; pathogenic/likely pathogenic and somatic-flagged filtered              |
| OncoKB          | Current      | API  | ~7,700                          | ~850 genes; 130 cancer types; academic token required                             |
| TP53 database   | R21          | File | ~29,900                         | GRCh38 CSV; functional annotations for >9,000 mutant proteins                     |
| CancerHotspots  | v2           | API  | ~3,181                          | Live REST API; 24,592 tumour samples; q-value filtered                             |
|                 |              |      |                                 |                                                                                    |
| **Raw total**   |              |      | **~33,391,744**                 | Pre-deduplication; significant inter-database overlap expected                     |

**Note on cBioPortal:** cBioPortal has been replaced by the TCGA mc3 PanCancer Atlas MAF as the
pan-cancer count source. The TCGA dataset provides broader, reproducible, offline coverage of the
same tumour population without live API dependency or 503 variability. `parse_cbioportal.py` is
retained for optional use but is disabled by default.

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
│   ├── parse_cbioportal.py  # retained; disabled by default
│   ├── parse_tcga.py
│   └── parse_hotspots.py
├── tools/                   # Standalone utilities; run manually
│   ├── annotate_panels.py
│   ├── clinvar_vep_annotate.py
│   ├── db_fix.py            # Liftover GRCh37->GRCh38 for GENIE and TCGA MAFs
│   ├── fetch_cbioportal.py  # Optional: snapshot cBioPortal API data
│   ├── hotspots_vep_remap.py
│   ├── mane_audit.py
│   └── mane_remap.py
└── logs/                    # Auto-created on first run
```

---

## Data preparation

GENIE and TCGA MAFs are distributed in GRCh37 and must be lifted to GRCh38
before the pipeline can use them. This is a one-time step per data version.

```bash
# 1. Download GENIE MAF from Synapse (requires free account and data agreement)
#    https://www.synapse.org/#!Synapse:syn7222066

# 2. Download TCGA mc3 MAF from GDC PanCancer Atlas
#    https://gdc.cancer.gov/about-data/publications/pancanatlas
#    File: mc3.v0.2.8.PUBLIC.maf.gz  (md5: 639ad8f8386e98dacc22e439188aa8fa)

# 3. Run liftover (lifts both GENIE and TCGA; writes logs alongside each file)
python3 tools/db_fix.py \
    --genie-maf data/genie/data_mutations_extended.txt \
    --tcga-maf  data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz
```

Liftover logs are written alongside each lifted file recording date, chain file,
row counts, and discard rate.

---

## Usage

### 1. Build the whitelist

Install Python dependencies first:
```bash
pip install -r requirements.txt
```

Then run the pipeline:
```bash
bash run_pipeline.sh /path/to/data/
```

Optional flags:
```bash
# Skip parsers; re-run merge/filter/output from saved intermediates
bash run_pipeline.sh --from-intermediates

# Run without OncoKB (no token required)
bash run_pipeline.sh --skip-sources oncokb

# Skip specific sources
bash run_pipeline.sh --skip-sources cosmic,tcga
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

HGVSc strings in the whitelist are carried through from source databases where
available. They are not re-annotated by ONCOSIEVE. Transcript concordance across
sources is not guaranteed; use coordinates (chrom, pos, ref, alt) as the
authoritative join key.

---

## Configuration

Edit `settings.yaml` to change thresholds, VAF floors, and OncoKB token.
Edit `config.yaml` to change file paths or disable sources.

---

## OncoKB token

Academic licence. Token expires in 6 months. Update in `settings.yaml` under `oncokb.api_token`.

Register at: https://www.oncokb.org/account/register

---

## References

If you use ONCOSIEVE in published work, please cite the underlying databases
using the references below.

**COSMIC**
Tate JG, Bamford S, Jubb HC, et al. COSMIC: the Catalogue Of Somatic Mutations
In Cancer. Nucleic Acids Res. 2019;47(D1):D941-D947.
doi:10.1093/nar/gky1015. PMID:30371878

**AACR Project GENIE**
The AACR Project GENIE Consortium. AACR Project GENIE: Powering Precision
Medicine through an International Consortium. Cancer Discov.
2017;7(8):818-831. doi:10.1158/2159-8290.CD-17-0151.
Include the version of GENIE data used (e.g. v19.0) in your methods.

**TCGA mc3 PanCancer Atlas**
Ellrott K, Bailey MH, Saksena G, et al. Scalable Open Science Approach for
Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines. Cell Syst.
2018;6(3):271-281.e7. doi:10.1016/j.cels.2018.03.002. PMID:29596782

**OncoKB**
Suehnholz SP, Nissan MH, Zhang H, et al. Quantifying the Expanding Landscape
of Clinical Actionability for Patients with Cancer. Cancer Discov.
2024;14(1):49-65. doi:10.1158/2159-8290.CD-23-0467. PMID:37849038

Chakravarty D, Gao J, Phillips SM, et al. OncoKB: A Precision Oncology
Knowledge Base. JCO Precis Oncol. 2017;1:PO.17.00011.
doi:10.1200/PO.17.00011. PMID:28890946

**ClinVar**
Landrum MJ, Lee JM, Benson M, et al. ClinVar: public archive of
interpretations of clinically relevant variants. Nucleic Acids Res.
2016;44(D1):D862-868. doi:10.1093/nar/gkv1222. PMID:26582918

**NCI TP53 Database**
de Andrade KC, Lee EE, Tookmanian EM, et al. The TP53 Database: transition
from the International Agency for Research on Cancer to the US National Cancer
Institute. Cell Death Differ. 2022;29(5):1071-1073.
doi:10.1038/s41418-022-00976-3. PMID:35352025

Database release: The TP53 Database (R21, Jan 2025): https://tp53.cancer.gov

**Cancer Hotspots**
Bandlamudi C, et al. Cancer type-specific variation in patterns of driver
alterations across 50,000 tumors. (2026). cancerhotspots.org

Chang MT, Asthana S, Gao SP, et al. Identifying recurrent mutations in cancer
reveals widespread lineage diversity and mutational specificity. Nat Biotechnol.
2016;34(2):155-163. doi:10.1038/nbt.3391. PMID:26619011

Chang MT, et al. Accelerating discovery of functional mutant alleles in cancer.
Cancer Discov. 2018;8(2):174-183. doi:10.1158/2159-8290.CD-17-0321.
PMID:29247016
