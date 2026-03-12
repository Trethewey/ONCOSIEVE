#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Mutect2 post-filter rescue: apply the ONCOSIEVE whitelist to rescued low-VAF variants
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
mutect2_rescue.py
Post-filter rescue of Mutect2-filtered variants using the pan-cancer whitelist.

This script takes:
  1. A Mutect2 VCF that has been processed through FilterMutectCalls
     (variants have FILTER values set, e.g. 'weak_evidence', 'low_allele_frac', etc.)
  2. The whitelist VCF produced by build_whitelist.py

For each variant in the Mutect2 VCF:
  - If FILTER = PASS: leave unchanged.
  - If FILTER != PASS (i.e. filtered-out): check whether the variant is
    in the whitelist AND the tumour VAF >= min_vaf_rescue (default 0.005).
    If both conditions are met: set FILTER = 'whitelist_rescued' and
    annotate with INFO/WL_TIER and INFO/WL_SOURCES.

Usage:
    python mutect2_rescue.py \\
        --input   sample.filtered.vcf.gz \\
        --whitelist pan_cancer_whitelist_GRCh38.vcf.gz \\
        --output  sample.rescued.vcf.gz \\
        --min-vaf 0.005 \\
        --tumor-sample TUMOR   # as named in the VCF FORMAT column

Output:
    A new VCF with rescued variants. Requires pysam >= 0.21.

Notes:
  - Variants rescued by the whitelist retain their original FORMAT fields
    (GT, AD, DP, AF etc.) unchanged. Only the FILTER column is modified.
  - The original FILTER value is preserved in INFO/ORIGINAL_FILTER.
  - Whitelist rescue is additive: variants already PASS are not modified.
  - Requires the input VCF to be bgzipped and tabix-indexed, OR uncompressed.
    A non-indexed gzip VCF is also handled (slower, no random access).
"""

import argparse
import gzip
import os
import sys
from typing import Optional

try:
    import pysam
except ImportError:
    sys.exit('ERROR: pysam is required.  Install with: pip install pysam')

import yaml

RESCUE_FILTER_LABEL = 'whitelist_rescued'


def load_whitelist(wl_vcf: str) -> dict[tuple, dict]:
    """
    Load whitelist VCF into an in-memory dict keyed by (chrom, pos, ref, alt).
    Value: dict with WL_TIER, WL_SOURCES, ONCOKB.

    Parameters
    ----------
    wl_vcf : path to whitelist VCF (plain, gzip, or bgzip+tabix)
    """
    wl: dict[tuple, dict] = {}
    vcf = pysam.VariantFile(wl_vcf, 'r')
    for rec in vcf.fetch():
        for alt in rec.alts or []:
            key = (rec.chrom, rec.pos, rec.ref, alt)
            wl[key] = {
                'wl_tier':    rec.info.get('WL_TIER', '.'),
                'wl_sources': rec.info.get('SOURCES', '.'),
                'oncokb':     rec.info.get('ONCOKB', ''),
                'n_samples':  rec.info.get('N_SAMPLES', 0),
            }
    vcf.close()
    print(f'[mutect2_rescue] Whitelist loaded: {len(wl):,} variants', flush=True)
    return wl


def get_tumour_sample_index(header: pysam.VariantHeader,
                             tumour_name: Optional[str]) -> int:
    """
    Return the index of the tumour sample in the VCF sample list.
    If tumour_name is not given, assumes the first sample is tumour.
    Raises ValueError if the name is not found.
    """
    samples = list(header.samples)
    if not samples:
        raise ValueError('VCF has no sample columns — cannot extract VAF.')
    if tumour_name is None:
        return 0
    if tumour_name not in samples:
        raise ValueError(
            f'Tumour sample "{tumour_name}" not found in VCF. '
            f'Available samples: {samples}'
        )
    return samples.index(tumour_name)


def get_vaf(record: pysam.VariantRecord, sample_idx: int) -> Optional[float]:
    """
    Extract tumour VAF from a Mutect2 record.
    Tries FORMAT/AF first, then computes from FORMAT/AD.
    Returns None if VAF cannot be determined.
    """
    sample = record.samples[sample_idx]

    # Mutect2 writes AF directly in FORMAT
    af = sample.get('AF')
    if af is not None:
        try:
            return float(af[0]) if hasattr(af, '__iter__') else float(af)
        except (TypeError, ValueError, IndexError):
            pass

    # Fallback: compute from AD
    ad = sample.get('AD')
    if ad is not None and len(ad) >= 2:
        ref_depth, alt_depth = ad[0], ad[1]
        total = ref_depth + alt_depth
        if total > 0:
            return alt_depth / total

    return None


def rescue(input_vcf: str,
           whitelist_vcf: str,
           output_vcf: str,
           min_vaf: float = 0.005,
           tumour_sample: Optional[str] = None,
           vaf_floors: dict | None = None) -> None:
    """
    Main rescue function.

    Parameters
    ----------
    input_vcf      : Mutect2 FilterMutectCalls output VCF (bgzip+tabix or plain)
    whitelist_vcf  : whitelist VCF from build_whitelist.py
    output_vcf     : output path (.vcf, .vcf.gz, or .vcf.bgz)
    min_vaf        : minimum VAF to rescue (inclusive)
    tumour_sample  : name of tumour sample column in the VCF; None = first sample
    """
    if vaf_floors is None:
        vaf_floors = {'tier1': min_vaf, 'tier2': min_vaf, 'tier3': min_vaf}
    wl = load_whitelist(whitelist_vcf)

    in_vcf  = pysam.VariantFile(input_vcf, 'r')
    header  = in_vcf.header.copy()

    # Add new FILTER and INFO entries to the header
    if RESCUE_FILTER_LABEL not in header.filters:
        header.add_line(
            f'##FILTER=<ID={RESCUE_FILTER_LABEL},'
            f'Description="Variant rescued by pan-cancer whitelist; '
            f'original filter overridden. VAF >= {min_vaf}.">'
        )
    if 'ORIGINAL_FILTER' not in header.info:
        header.add_line(
            '##INFO=<ID=ORIGINAL_FILTER,Number=.,Type=String,'
            'Description="Original Mutect2 FILTER value before whitelist rescue">'
        )
    if 'WL_TIER' not in header.info:
        header.add_line(
            '##INFO=<ID=WL_TIER,Number=1,Type=Integer,'
            'Description="Whitelist tier (1=highest confidence, 3=minimum threshold)">'
        )
    if 'WL_SOURCES' not in header.info:
        header.add_line(
            '##INFO=<ID=WL_SOURCES,Number=1,Type=String,'
            'Description="Whitelist sources supporting this variant">'
        )
    if 'WL_N_SAMPLES' not in header.info:
        header.add_line(
            '##INFO=<ID=WL_N_SAMPLES,Number=1,Type=Integer,'
            'Description="Number of samples in whitelist databases">'
        )
    if 'WL_ONCOKB' not in header.info:
        header.add_line(
            '##INFO=<ID=WL_ONCOKB,Number=1,Type=String,'
            'Description="OncoKB oncogenicity classification from whitelist">'
        )

    tumour_idx = get_tumour_sample_index(header, tumour_sample)
    print(f'[mutect2_rescue] Tumour sample index: {tumour_idx} '
          f'({list(header.samples)[tumour_idx]})', flush=True)

    out_mode = 'wz' if output_vcf.endswith('.gz') or \
                        output_vcf.endswith('.bgz') else 'w'
    out_vcf  = pysam.VariantFile(output_vcf, out_mode, header=header)

    n_total = n_already_pass = n_eligible = n_rescued = n_vaf_fail = 0

    for rec in in_vcf.fetch():
        n_total += 1
        filters = list(rec.filter.keys())

        # Already passing — write unchanged
        if filters == ['PASS'] or filters == []:
            n_already_pass += 1
            out_vcf.write(rec)
            continue

        # Variant is filtered: check whitelist
        for alt in rec.alts or []:
            key = (rec.chrom, rec.pos, rec.ref, alt)
            wl_entry = wl.get(key)

            if wl_entry is None:
                continue   # not in whitelist

            n_eligible += 1

            # Tier-aware VAF floor
            tier = int(wl_entry.get('wl_tier', 3)) if str(wl_entry.get('wl_tier', 3)).isdigit() else 3
            if tier == 1:
                vaf_floor = vaf_floors.get('tier1', min_vaf)
            elif tier == 2:
                vaf_floor = vaf_floors.get('tier2', min_vaf)
            else:
                vaf_floor = vaf_floors.get('tier3', min_vaf)

            vaf = get_vaf(rec, tumour_idx)
            if vaf is None or vaf < vaf_floor:
                n_vaf_fail += 1
                break

            # RESCUE: modify the record
            original_filters = '|'.join(filters)
            rec.filter.clear()
            rec.filter.add(RESCUE_FILTER_LABEL)
            rec.info['ORIGINAL_FILTER'] = original_filters
            rec.info['WL_TIER']         = int(wl_entry['wl_tier']) \
                                          if str(wl_entry['wl_tier']).isdigit() else 3
            rec.info['WL_SOURCES']      = str(wl_entry['wl_sources'])
            rec.info['WL_N_SAMPLES']    = int(wl_entry['n_samples'])
            if wl_entry.get('oncokb'):
                rec.info['WL_ONCOKB']   = str(wl_entry['oncokb'])

            n_rescued += 1
            break   # one matching alt is enough to rescue the record

        out_vcf.write(rec)

    in_vcf.close()
    out_vcf.close()

    print(f'\n[mutect2_rescue] Summary', flush=True)
    print(f'  Total input variants : {n_total:,}')
    print(f'  Already PASS         : {n_already_pass:,}')
    print(f'  In whitelist         : {n_eligible:,}')
    print(f'  VAF < {min_vaf:.1%} (not rescued) : {n_vaf_fail:,}')
    print(f'  Rescued              : {n_rescued:,}')
    print(f'  Output written to    : {output_vcf}')

    if output_vcf.endswith('.gz') or output_vcf.endswith('.bgz'):
        print(f'\n  Index with: tabix -p vcf {output_vcf}')


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Rescue Mutect2-filtered variants using a pan-cancer whitelist.'
    )
    ap.add_argument('--input',    required=True,
                    help='Mutect2 filtered VCF (FilterMutectCalls output)')
    ap.add_argument('--whitelist', required=True,
                    help='Whitelist VCF from build_whitelist.py')
    ap.add_argument('--output',   required=True,
                    help='Output VCF path (.vcf or .vcf.gz)')
    ap.add_argument('--min-vaf',  type=float, default=0.005,
                    help='Minimum tumour VAF to rescue (default: 0.005 = 0.5%%)')
    ap.add_argument('--tumor-sample', default=None,
                    help='Name of the tumour sample column in the VCF '
                         '(default: first sample)')
    ap.add_argument('--config', default=None,
                    help='Optional config.yaml to read min_vaf_rescue from '
                         '(overrides --min-vaf)')
    args = ap.parse_args()

    min_vaf = args.min_vaf
    if args.config and os.path.exists(args.config):
        with open(args.config) as fh:
            cfg = yaml.safe_load(fh)
        min_vaf = cfg.get('thresholds', {}).get('min_vaf_rescue', min_vaf)
        settings_file = cfg.get('settings_file', 'settings.yaml')
    else:
        cfg = {}
        settings_file = 'settings.yaml'

    # Load tier-aware VAF floors from settings.yaml if available
    vaf_floors = {'tier1': min_vaf, 'tier2': min_vaf, 'tier3': min_vaf}
    if os.path.exists(settings_file):
        with open(settings_file) as fh:
            settings = yaml.safe_load(fh) or {}
        vaf_floors = settings.get('vaf_rescue', vaf_floors)

    print(f'[mutect2_rescue] VAF rescue floors: '
          f'Tier1={vaf_floors.get("tier1",min_vaf):.3f}  '
          f'Tier2={vaf_floors.get("tier2",min_vaf):.3f}  '
          f'Tier3={vaf_floors.get("tier3",min_vaf):.3f}')

    rescue(
        input_vcf      = args.input,
        whitelist_vcf  = args.whitelist,
        output_vcf     = args.output,
        min_vaf        = min_vaf,
        tumour_sample  = args.tumor_sample,
        vaf_floors     = vaf_floors,
    )


if __name__ == '__main__':
    main()
