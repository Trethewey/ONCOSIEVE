# Changelog

All notable changes to OncoSieve are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] — 2026-05-19

### Added
- **PrimateAI-3D annotation** (`tools/annotate_primateai.py`). Polars lazy-scan over the 1.7 GB hg38 dataset; adds `primateai_score`, `primateai_percentile`, and `primateai_prediction` columns to the post-pipeline output. Pipeline auto-detects `data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz`.
- **PrimateAI-3D citation** in the references table (Gao et al., *Science* 2023).
- **Sources table** in the HTML report, auto-generated from `output/database_versions.txt`. Shows database, version/date, file vs API, and origin.
- **CSV download** button on the high-confidence variant table (DataTables Buttons).
- **Favicon** in the HTML report, embedded as a base64 data URI for self-contained portability.
- **Packaging directory** with `Dockerfile`, `docker-compose.yml`, conda `environment.yml`, `setup.py` (pip), and `Makefile`.
- **MIT licence** (`LICENSE`); previously labelled "Research Use Only" in `setup.py`.
- **CHANGELOG.md** (this file).

### Changed
- **Brand:** `ONCOSIEVE` → `OncoSieve` across code surfaces, HTML report, Excel sheet titles, and `setup.py`. File paths and the GitHub repo URL retain their casing.
- **Repository renamed** on GitHub from `ONCOSIEVE` to `OncoSieve`. Old URL continues to redirect.
- **HTML report header redesigned:** mark SVG instead of cropped raster; transparent background; new wordmark in `-apple-system/Segoe UI` matching the SVG logo; "2.0" version badge baseline-aligned beside the title; restructured subtitle (`Pan-Cancer Variant Whitelist – Generated DATE`); source line lists databases + REVEL + PrimateAI-3D.
- **HTML report variant table** now uses DataTables `data:` JSON-init instead of pre-rendered `<tr>` rows. ~40% smaller HTML (6.8 → 4.1 MB) and renders effectively instantly.
- **High-Confidence KPI card** recoloured from teal-green to logo red `#C0392B`.
- **References section** in the report: PrimateAI-3D entry added; text colour grey → white; font size reduced by 2 pt.
- **TP53 database version** pinned in `config.yaml` (`tp53.version: R21`); `build_whitelist.py` reads this in preference to the previous mtime fallback.
- **`run_oncosieve.sh`** auto-detects `data/PRIMATE_AI/PrimateAI-3D.hg38.txt.gz` and passes `--primateai` to `post_pipeline.py` when present.
- **`.gitignore`** reorganised into labelled sections; adds `settings.yaml`, `NOTES.md`, `nohup.out`.
- **README** rewritten: new MIT badge under the logo, Annotation Layers sub-table, updated source versions (ClinVar 2026-03-09, OncoKB v7.1 Apr 2026), new PrimateAI-3D data-prep section, new Packaging section, Licence section, output columns 23–25 documented, tiering section accurately states PrimateAI-3D is annotation-only in v2.
- **`settings.yaml`** untracked. The template is `settings.yaml.example` (tracked); users copy it to `settings.yaml` and add their OncoKB token locally.
- **VCF source header** simplified to `pan_cancer_whitelist_pipeline` (no version suffix).

### Removed
- **`packaging/install_and_run.bat`** — the Windows .bat installer already depended on Git Bash / WSL (for `tee` and bcftools). Docker is now the supported Windows distribution path.
- **Tool-version strings** (`ONCOSIEVE v1.0`, `pipeline_v1.0`) from code surfaces, banner, Excel title, and VCF source header. The version lives in the Git tag and the in-report "2.0" badge.
- **Tracked `settings.yaml`** — see Changed.

### Not implemented in this release
- **PrimateAI-3D tier promotion.** Documented as a future enhancement in `tools/annotate_primateai.py:326`. The scores are annotated but do not feed into the tier assignment.

## [1.0.0] — 2026-03-16

### Added
- Initial release of OncoSieve.
- Pan-cancer somatic variant whitelist built from seven curated databases:
  COSMIC v103, AACR GENIE v19.0, TCGA mc3 v0.2.8, ClinVar, OncoKB,
  TP53 Database, CancerHotspots v2 — ~46.4 million raw variants pre-deduplication.
- Tiered output (Tier 1 / 2 / 3) with VAF floors configurable per tier.
- MANE Select enrichment for transcript IDs.
- Mutect2 rescue (`mutect2_rescue.py`): applies the whitelist to FilterMutectCalls
  VCFs with tier-aware VAF floors.
- Post-pipeline (`tools/post_pipeline.py`): REVEL annotation, high-confidence
  filtering (drops ClinVar-only zero-sample rows), Excel exports.
- Interactive HTML report with Plotly charts and a DataTables variant browser.
- Orchestrator `run_oncosieve.sh` with `--from-intermediates`, `--skip-sources`,
  `--rescue` flags.
- Pre-run audit (`pre_check.py`) and validation harness (`test_pipeline.py`).
