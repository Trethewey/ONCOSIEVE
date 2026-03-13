#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — post_process_whitelist.py
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
#
# Run after build_whitelist.py and mane_remap.py.
#
# Adds three columns to the final whitelist TSV:
#
#   genome_version  — The genome build for each row's coordinate set.
#                     Derived from the sources column. All sources should be
#                     GRCh38 after db_fix.py liftover. Any residual GRCh37
#                     rows are flagged (pass --sources-grch37 to name them).
#
#   transcript_id   — ENST accession extracted from hgvsc, version stripped.
#                     E.g. ENST00000288602 from ENST00000288602.6:c.1799T>A.
#                     RefSeq NM_ accessions are also extracted.
#                     Set to '' if hgvsc is absent or unparseable.
#
#   protein_change  — Normalised protein change in 3-letter HGVS amino acid
#                     format (e.g. p.Val600Glu). Handles:
#                       1-letter input:  p.V600E  -> p.Val600Glu
#                       Lowercase:       p.val600glu -> p.Val600Glu
#                       No prefix:       V600E -> p.Val600Glu
#                       Stop-gained:     p.R213* -> p.Arg213Ter
#                       Frameshift:      p.Ala123fs, p.Ala123fsTer45
#                       In-frame indel:  p.Glu746_Ala750del, p.T599_V600insE
#                       Synonymous:      p.= or p.Val600=
#                     Set to '' if hgvsp is absent or unparseable.
#
# NOTE on HGVSc/HGVSp and liftover:
#   These fields are transcript-relative and are NOT modified by coordinate
#   liftover. protein_change is derived from hgvsp, which remains valid
#   across builds. If ref allele mismatches were flagged by db_fix.py,
#   re-annotate those rows with VEP before running this script.
#
# Usage:
#   python3 post_process_whitelist.py \
#       --whitelist pan_cancer_whitelist_GRCh38.tsv.gz \
#       --out       pan_cancer_whitelist_GRCh38_annotated.tsv.gz
#
# Optional:
#   --sources-grch37  Comma-separated source names still carrying GRCh37
#                     coordinates (post-liftover failure). Default: none.
#   --strict-hgvsp    Exit non-zero if any hgvsp could not be normalised.
# =============================================================================

import argparse
import gzip
import re
import sys
from typing import Optional


# ── Amino acid tables ─────────────────────────────────────────────────────────

AA1_TO_AA3 = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    'X': 'Ter', '*': 'Ter', '?': '?',
}

AA3_CANON = {v.upper(): v for v in AA1_TO_AA3.values()}
AA3_SET   = {v.upper() for v in AA1_TO_AA3.values()}


def _aa_to_3(aa: str) -> str:
    aa = aa.strip()
    if len(aa) == 1:
        return AA1_TO_AA3.get(aa.upper(), aa)
    if len(aa) == 3 and aa.upper() in AA3_SET:
        return AA3_CANON.get(aa.upper(), aa.capitalize())
    return aa


# ── Regex patterns ────────────────────────────────────────────────────────────

_RE_MISSENSE   = re.compile(r'^([A-Za-z*]{1,3})(\d+)([A-Za-z*]{1,3})$')
_RE_NONSENSE   = re.compile(r'^([A-Za-z*]{1,3})(\d+)(?:\*|Ter|ter|TER|X)$')
_RE_FRAMESHIFT = re.compile(r'^([A-Za-z*]{1,3})(\d+)fs(?:\*|Ter|ter)?(\d*)$', re.I)
_RE_DEL        = re.compile(r'^([A-Za-z*]{1,3})(\d+)(?:_([A-Za-z*]{1,3})(\d+))?del$', re.I)
_RE_INS        = re.compile(r'^([A-Za-z*]{1,3})(\d+)_([A-Za-z*]{1,3})(\d+)ins([A-Za-z*]+)$', re.I)
_RE_DUP        = re.compile(r'^([A-Za-z*]{1,3})(\d+)(?:_([A-Za-z*]{1,3})(\d+))?dup$', re.I)
_RE_SYN        = re.compile(r'^(?:([A-Za-z*]{1,3})(\d+))?=$')
_RE_EXT        = re.compile(r'^\*(\d+)([A-Za-z*]{1,3})ext(?:\*|Ter)?(\d*)$', re.I)
_RE_ENST       = re.compile(r'(ENST\d+)', re.I)


# ── Protein change normalisation ──────────────────────────────────────────────

def normalise_protein_change(raw: str) -> str:
    """
    Normalise a protein change string to 3-letter HGVS format.
    Returns '' for absent or unparseable input.
    """
    if not raw or not isinstance(raw, str):
        return ''
    raw = raw.strip()
    if not raw or raw in {'.', '-', 'NA', 'N/A', 'nan', 'None'}:
        return ''

    # Strip ENSP prefix e.g. ENSP00000484643.1:p.Trp235Arg -> Trp235Arg
    raw = re.sub(r'^ENSP\d+(?:\.\d+)?:p\.', '', raw, flags=re.I)
    # Strip plain ENSP prefix without p. e.g. ENSP00000484643.1:Trp235Arg
    raw = re.sub(r'^ENSP\d+(?:\.\d+)?:', '', raw, flags=re.I)

    body = raw[2:] if raw.lower().startswith('p.') else raw

    if body == '?':
        return 'p.?'

    m = _RE_SYN.match(body)
    if m:
        if m.group(1):
            return f'p.{_aa_to_3(m.group(1))}{m.group(2)}='
        return 'p.='

    # Nonsense before missense — patterns overlap at the stop character
    m = _RE_NONSENSE.match(body)
    if m:
        return f'p.{_aa_to_3(m.group(1))}{m.group(2)}Ter'

    m = _RE_FRAMESHIFT.match(body)
    if m:
        aa, pos, ext = m.group(1), m.group(2), m.group(3)
        return f'p.{_aa_to_3(aa)}{pos}fs{"Ter" + ext if ext else ""}'

    m = _RE_DEL.match(body)
    if m:
        aa1, pos1 = _aa_to_3(m.group(1)), m.group(2)
        if m.group(3):
            return f'p.{aa1}{pos1}_{_aa_to_3(m.group(3))}{m.group(4)}del'
        return f'p.{aa1}{pos1}del'

    m = _RE_DUP.match(body)
    if m:
        aa1, pos1 = _aa_to_3(m.group(1)), m.group(2)
        if m.group(3):
            return f'p.{aa1}{pos1}_{_aa_to_3(m.group(3))}{m.group(4)}dup'
        return f'p.{aa1}{pos1}dup'

    m = _RE_INS.match(body)
    if m:
        aa1, pos1 = _aa_to_3(m.group(1)), m.group(2)
        aa2, pos2 = _aa_to_3(m.group(3)), m.group(4)
        ins_aas   = ''.join(_aa_to_3(c) for c in
                            re.findall(r'[A-Za-z*]{1,3}', m.group(5)))
        return f'p.{aa1}{pos1}_{aa2}{pos2}ins{ins_aas}'

    m = _RE_EXT.match(body)
    if m:
        pos, aa, ext = m.group(1), _aa_to_3(m.group(2)), m.group(3)
        return f'p.*{pos}{aa}ext{"Ter" + ext if ext else ""}'

    m = _RE_MISSENSE.match(body)
    if m:
        return f'p.{_aa_to_3(m.group(1))}{m.group(2)}{_aa_to_3(m.group(3))}'

    # Return raw with p. prefix if nothing matched
    return f'p.{body}' if not raw.lower().startswith('p.') else raw


# ── Transcript ID extraction ──────────────────────────────────────────────────

def extract_transcript_id(hgvsc: str) -> str:
    if not hgvsc or not isinstance(hgvsc, str) or hgvsc in {'.', '-'}:
        return ''
    m = _RE_ENST.search(hgvsc)
    if m:
        return m.group(1).upper()
    m = re.match(r'^(N[MRP]_\d+)', hgvsc, re.I)
    if m:
        return m.group(1).upper()
    return ''


# ── Genome version assignment ─────────────────────────────────────────────────

def assign_genome_version(sources: str, grch37_sources: set) -> str:
    if not grch37_sources or not sources or sources in {'.', '-', ''}:
        return 'GRCh38'
    src_set = {s.strip().lower() for s in re.split(r'[,;|]', sources) if s.strip()}
    if src_set and src_set.issubset({s.lower() for s in grch37_sources}):
        return 'GRCh37'
    return 'GRCh38'


# ── I/O helpers ───────────────────────────────────────────────────────────────

def _open_r(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'rt')

def _open_w(path):
    return gzip.open(path, 'wt') if path.endswith('.gz') else open(path, 'wt')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='ONCOSIEVE post_process_whitelist.py'
    )
    parser.add_argument('--whitelist',      required=True)
    parser.add_argument('--out',            required=True)
    parser.add_argument('--sources-grch37', default='')
    parser.add_argument('--strict-hgvsp',   action='store_true')
    args = parser.parse_args()

    grch37_set: set = set()
    if args.sources_grch37:
        grch37_set = {s.strip() for s in args.sources_grch37.split(',') if s.strip()}

    print(f'Input:   {args.whitelist}')
    print(f'Output:  {args.out}')
    if grch37_set:
        print(f'GRCh37 sources flagged: {grch37_set}')
    print()

    new_cols = ['genome_version', 'transcript_id', 'protein_change']
    n_total = n_grch37 = n_no_transcript = n_no_protein = n_unparsed = 0
    unparsed_examples: list = []

    with _open_r(args.whitelist) as fh_in, _open_w(args.out) as fh_out:
        header_written = False
        col_i: dict = {}
        for line in fh_in:
            row = line.rstrip('\n').split('\t')
            if not header_written:
                col_i = {c: i for i, c in enumerate(row)}
                missing = [c for c in ['sources', 'hgvsc', 'hgvsp']
                           if c not in col_i]
                if missing:
                    print(f'ERROR: Missing required columns: {missing}',
                          file=sys.stderr)
                    sys.exit(1)
                fh_out.write('\t'.join(row + new_cols) + '\n')
                header_written = True
                continue

            sources = row[col_i.get('sources', -1)] if 'sources' in col_i else ''
            hgvsc   = row[col_i.get('hgvsc',   -1)] if 'hgvsc'   in col_i else ''
            hgvsp   = row[col_i.get('hgvsp',   -1)] if 'hgvsp'   in col_i else ''

            genome_version = assign_genome_version(sources, grch37_set)
            transcript_id  = extract_transcript_id(hgvsc)
            protein_change = normalise_protein_change(hgvsp)

            n_total += 1
            if genome_version == 'GRCh37':
                n_grch37 += 1
            if not transcript_id:
                n_no_transcript += 1
            if not protein_change:
                n_no_protein += 1
            # Flag if input hgvsp was present but normalisation left a
            # non-standard string (no 'p.' prefix after normalisation)
            if (hgvsp and hgvsp not in {'.', '-', 'NA', 'nan', ''}
                    and protein_change
                    and not protein_change.startswith('p.')):
                n_unparsed += 1
                if len(unparsed_examples) < 10:
                    unparsed_examples.append(hgvsp)

            fh_out.write('\t'.join(row + [genome_version, transcript_id,
                                          protein_change]) + '\n')

    print(f'Rows processed:         {n_total:,}')
    print(f'  GRCh38:               {n_total - n_grch37:,}')
    if n_grch37:
        print(f'  GRCh37 flagged:       {n_grch37:,}')
    print(f'  With transcript ID:   {n_total - n_no_transcript:,}')
    print(f'  Without transcript:   {n_no_transcript:,}  '
          '(OncoKB/ClinVar concept-level variants expected here)')
    print(f'  With protein change:  {n_total - n_no_protein:,}')
    print(f'  No protein change:    {n_no_protein:,}  '
          '(splice/intronic/noncoding rows expected)')
    if n_unparsed:
        print(f'  Unparsed hgvsp:       {n_unparsed:,}')
        print('  Examples (first 10):')
        for ex in unparsed_examples:
            print(f'    {ex}')
        if args.strict_hgvsp:
            sys.exit(1)

    print(f'\nOutput written: {args.out}')


if __name__ == '__main__':
    main()
