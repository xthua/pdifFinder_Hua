#!/usr/bin/env python3
"""
pdif_pwm_integration.py - Integration layer for PWM-based pdif pwm scanner.

Provides clean entry points for loading PWMs, scanning sequences,
deduplicating hits, and comparing with the old pdifFinder.

Author: PdifFinder Team
License: GPLv3
"""

import os
import sys
import json
import subprocess
import tempfile
import shutil
from collections import OrderedDict

from Bio import SeqIO
from Bio.Seq import Seq

from PdifFinder.pdif_pwm_scanner import (
    build_pwms, load_pwms, scan_sequence,
    ARM_LENGTH, SPACER_LENGTH, SITE_LENGTH
)

# ---------------------------------------------------------------------------
# Module-level PWM cache
# ---------------------------------------------------------------------------

_pwm_cache = None


def _get_pwms():
    """Return cached (xerC_pwm, xerD_pwm, xerD_degen_pwm, xerC_seeds, xerD_seeds)."""
    global _pwm_cache
    if _pwm_cache is None:
        json_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'data', 'xer_pwm.json'
        )
        if not os.path.exists(json_path):
            build_pwms()
        _pwm_cache = load_pwms(json_path)
    return _pwm_cache


def clear_pwm_cache():
    """Clear the cached PWMs (useful for testing)."""
    global _pwm_cache
    _pwm_cache = None


# ---------------------------------------------------------------------------
# Sequence loading
# ---------------------------------------------------------------------------

def load_fasta_sequence(filepath):
    """Load the first sequence from a FASTA file.

    Returns
    -------
    tuple of (seq_id, sequence_string) or (None, None)
    """
    with open(filepath, 'r') as fh:
        for rec in SeqIO.parse(fh, 'fasta'):
            return rec.id, str(rec.seq).upper()
    return None, None


def load_genbank_sequence(filepath):
    """Load the first sequence from a GenBank file.

    Returns
    -------
    tuple of (seq_id, sequence_string) or (None, None)
    """
    with open(filepath, 'r') as fh:
        for rec in SeqIO.parse(fh, 'genbank'):
            return rec.id, str(rec.seq).upper()
    return None, None


def load_sequence(filepath):
    """Auto-detect format and load first sequence.

    Returns
    -------
    tuple of (seq_id, sequence_string) or (None, None)
    """
    base = os.path.splitext(filepath)[1].lower()
    if base in ('.gb', '.genbank', '.gbk'):
        return load_genbank_sequence(filepath)
    return load_fasta_sequence(filepath)


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def deduplicate_sites(sites, bp=2):
    """Merge overlapping or nearby pdif sites within *bp* bases.

    When two hits overlap or are within *bp* positions, the one with the
    higher total_score is kept.

    Parameters
    ----------
    sites : list of dict
        Results from ``scan_sequence`` (must have 'start' and 'total_score').
    bp : int
        Merge distance in bases (default 2).

    Returns
    -------
    list of dict
    """
    if not sites:
        return []
    sorted_sites = sorted(sites, key=lambda s: s['start'])
    merged = []
    for site in sorted_sites:
        if merged and site['start'] - merged[-1]['start'] <= bp:
            if site.get('total_score', 0.0) > merged[-1].get('total_score', 0.0):
                merged[-1] = site
        else:
            merged.append(site)
    return merged


# ---------------------------------------------------------------------------
# Scanning entry points
# ---------------------------------------------------------------------------

def scan_file(filepath, dedup_bp=2):
    """Scan a FASTA or GenBank file for pdif sites using pwm scanner.

    Parameters
    ----------
    filepath : str
        Path to FASTA or GenBank file.
    dedup_bp : int
        Merge distance for deduplication.

    Returns
    -------
    tuple of (seq_id, sites_list)
        sites_list is list of dict from scan_sequence, deduplicated.
    """
    seq_id, seq = load_sequence(filepath)
    if seq is None:
        return None, []
    xerC_pwm, xerD_pwm, xerD_degen, xerC_seeds, xerD_seeds = _get_pwms()
    raw_sites = scan_sequence(seq, xerC_pwm, xerD_pwm, xerC_seeds, xerD_seeds)
    sites = deduplicate_sites(raw_sites, bp=dedup_bp)
    return seq_id, sites


def scan_string(seq_str, seq_id='seq', dedup_bp=2):
    """Scan a DNA string for pdif sites.

    Parameters
    ----------
    seq_str : str
        DNA sequence.
    seq_id : str
        Identifier for the sequence.
    dedup_bp : int
        Merge distance for deduplication.

    Returns
    -------
    tuple of (seq_id, sites_list)
    """
    xerC_pwm, xerD_pwm, xerD_degen, xerC_seeds, xerD_seeds = _get_pwms()
    raw_sites = scan_sequence(seq_str, xerC_pwm, xerD_pwm, xerC_seeds, xerD_seeds)
    sites = deduplicate_sites(raw_sites, bp=dedup_bp)
    return seq_id, sites


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_sites(sites):
    """Classify pdif sites into CD and DC based on raw_orientation.

    Parameters
    ----------
    sites : list of dict

    Returns
    -------
    dict with keys: total, cd_sites, dc_sites, dc_proportion
    """
    total = len(sites)
    cd_sites = [s for s in sites if s.get('raw_orientation') == 'CD']
    dc_sites = [s for s in sites if s.get('raw_orientation') == 'DC']
    cd_count = len(cd_sites)
    dc_count = len(dc_sites)
    dc_prop = dc_count / total if total > 0 else 0.0
    return {
        'total': total,
        'cd': cd_count,
        'dc': dc_count,
        'dc_proportion': dc_prop,
        'cd_sites': cd_sites,
        'dc_sites': dc_sites,
    }


# ---------------------------------------------------------------------------
# Old pdifFinder runner
# ---------------------------------------------------------------------------

def run_old_pdif_finder(fasta_path, output_dir=None):
    """Run the legacy pdifFinder on a FASTA file (pdif_only mode, mocked AMR).

    Uses ``unittest.mock.patch`` to bypass the AMR gene guard so pdif
    detection runs without real BLAST hits.

    Parameters
    ----------
    fasta_path : str
        Path to input FASTA.
    output_dir : str or None
        Output directory (temp dir created if None).

    Returns
    -------
    (output_dir_path, return_code)
    """
    from unittest.mock import patch

    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='old_pdif_')
    else:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    import PdifFinder.pdifFinder as pf
    scripts_dir = os.path.dirname(os.path.abspath(pf.__file__))

    # Prepare directory structure
    os.makedirs(output_dir, exist_ok=True)
    tmp_dir = os.path.join(output_dir, 'tmp')
    os.makedirs(tmp_dir, exist_ok=True)
    for sub in ('pairSearch', 'seedSearch'):
        os.makedirs(os.path.join(tmp_dir, sub), exist_ok=True)

    # Placeholder files required by downstream steps
    for placeholder in ('AMRgene1.txt', 'AMRgene.txt'):
        ppath = os.path.join(output_dir, placeholder)
        if not os.path.exists(ppath):
            open(ppath, 'w').close()

    mock_find_rg = patch(
        'PdifFinder.pdifFinder.findResistanceGene',
        return_value=['1-1000'],
    )

    try:
        mock_find_rg.start()
        pf.singleThread(
            inFile=fasta_path,
            outdir=output_dir,
            blastnPath='mock-blastn',
            scriptsDir=scripts_dir,
            circle='true',
        )
        rc = 0
    except Exception as e:
        rc = -1
    finally:
        mock_find_rg.stop()

    return output_dir, rc


def _parse_pdif_site1(output_dir):
    """Parse pdif_site1.txt (raw output before finalFilter).

    Format varies; columns are: pdifName, start, end, xerC, spacer, xerD, [orientation, ...]
    """
    sites = []
    filepath = os.path.join(output_dir, 'pdif_site1.txt')
    if not os.path.exists(filepath):
        return sites
    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            orient = 'CD'
            if len(parts) >= 7:
                orient_raw = parts[6].strip()
                if orient_raw in ('C|D', 'CD'):
                    orient = 'CD'
                elif orient_raw in ('D|C', 'DC'):
                    orient = 'DC'
            sites.append({
                'name': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'xerC': parts[3],
                'spacer': parts[4],
                'xerD': parts[5],
                'orientation': orient,
                'pdif_db_name': '',
            })
    return sites


def parse_old_pdif_output(output_dir):
    """Parse the old pdifFinder output (pdif_site.txt or pdif_site1.txt).

    Format: seq_id | pdifName | start | end | xerC | spacer | xerD | orientation | pdif_db_name

    Returns
    -------
    list of dict with keys: name, start, end, xerC, spacer, xerD, orientation
    """
    pdif_site_file = os.path.join(output_dir, 'pdif_site.txt')
    if os.path.exists(pdif_site_file):
        sites = []
        with open(pdif_site_file, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 8:
                    orientation_raw = parts[7]
                    if orientation_raw == 'C|D':
                        orientation = 'CD'
                    elif orientation_raw == 'D|C':
                        orientation = 'DC'
                    else:
                        orientation = orientation_raw
                    sites.append({
                        'name': parts[1],
                        'start': int(parts[2]),
                        'end': int(parts[3]),
                        'xerC': parts[4],
                        'spacer': parts[5],
                        'xerD': parts[6],
                        'orientation': orientation,
                        'pdif_db_name': parts[8] if len(parts) >= 9 else '',
                    })
        return sites
    return _parse_pdif_site1(output_dir)


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

def compare_plasmid(plasmid_name, fasta_path,
                    output_dir=None,
                    download_cache_dir=None):
    """Run both old and new scanners on a plasmid, returning comparison stats.

    Parameters
    ----------
    plasmid_name : str
    fasta_path : str
    output_dir : str or None
    download_cache_dir : str or None - cache dir for downloads

    Returns
    -------
    dict with keys: name, old_total, old_cd, old_dc, old_dc_pct,
                    new_total, new_cd, new_dc, new_dc_pct
    """
    result = {
        'name': plasmid_name,
        'old_total': 0, 'old_cd': 0, 'old_dc': 0, 'old_dc_pct': 0.0,
        'new_total': 0, 'new_cd': 0, 'new_dc': 0, 'new_dc_pct': 0.0,
        'error': None,
    }

    # --- Old scanner ---
    try:
        old_outdir, rc = run_old_pdif_finder(fasta_path, output_dir)
        old_sites = parse_old_pdif_output(old_outdir)
        result['old_total'] = len(old_sites)
        result['old_cd'] = sum(1 for s in old_sites if s['orientation'] == 'CD')
        result['old_dc'] = sum(1 for s in old_sites if s['orientation'] == 'DC')
        if result['old_total'] > 0:
            result['old_dc_pct'] = round(result['old_dc'] / result['old_total'] * 100, 1)
        # Clean up old output dir
        if output_dir is None:
            shutil.rmtree(old_outdir, ignore_errors=True)
    except Exception as e:
        result['error'] = 'old: ' + str(e)

    # --- New pwm scanner ---
    try:
        seq_id, new_sites = scan_file(fasta_path)
        classification = classify_sites(new_sites)
        result['new_total'] = classification['total']
        result['new_cd'] = classification['cd']
        result['new_dc'] = classification['dc']
        if classification['total'] > 0:
            result['new_dc_pct'] = round(classification['dc_proportion'] * 100, 1)
    except Exception as e:
        if result['error']:
            result['error'] += '; new: ' + str(e)
        else:
            result['error'] = 'new: ' + str(e)

    return result


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

NCBI_EFETCH_URL = (
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    '?db=nucleotide&id={accession}&rettype=fasta&retmode=text'
)

CACHE_DIR = None


def _get_cache_dir():
    global CACHE_DIR
    if CACHE_DIR is None:
        CACHE_DIR = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            '..', 'tests', 'benchmark_data', 'plasmids'
        )
    return CACHE_DIR


def download_plasmid(accession, cache_dir=None, force=False):
    """Download a plasmid from NCBI if not cached.

    Parameters
    ----------
    accession : str
    cache_dir : str or None
    force : bool

    Returns
    -------
    str - path to cached FASTA file, or None on failure
    """
    if cache_dir is None:
        cache_dir = _get_cache_dir()
    elif not os.path.isdir(cache_dir):
        cache_dir = _get_cache_dir()

    # Check cache first
    fname = accession.replace('.', '_') + '.fasta'
    fname_alt = accession.replace('.', '_') + '_1.fasta'
    for cand in [fname, fname_alt]:
        path = os.path.join(cache_dir, cand)
        if os.path.exists(path):
            # Verify it's a FASTA file
            with open(path, 'r') as fh:
                header = fh.read(200).strip()
                if header.startswith('>'):
                    return path

    # Not cached - download
    url = NCBI_EFETCH_URL.format(accession=accession)
    outpath = os.path.join(cache_dir, fname)
    cmd = [
        'curl', '-s', '-k', '-L',
        '--connect-timeout', '30',
        '--max-time', '120',
        '-o', outpath,
        url
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=130)
        if os.path.exists(outpath) and os.path.getsize(outpath) > 0:
            # Check header
            with open(outpath, 'r') as fh:
                header = fh.read(200).strip()
                if header.startswith('>'):
                    return outpath
            # Not valid - remove
            os.remove(outpath)
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# DC proportion analysis (across all cached plasmids)
# ---------------------------------------------------------------------------

def analyze_all_cached_plasmids(cache_dir=None):
    """Run pwm scanner on ALL cached plasmid FASTA files and report DC stats.

    Returns
    -------
    list of dict with keys: name, total, cd, dc, dc_proportion
    """
    if cache_dir is None:
        cache_dir = _get_cache_dir()
    results = []
    for fname in sorted(os.listdir(cache_dir)):
        fpath = os.path.join(cache_dir, fname)
        if not os.path.isfile(fpath):
            continue
        if not (fname.endswith('.fasta') or fname.endswith('.fa')):
            continue
        # Skip all_plasmids.fasta (aggregated file)
        if fname == 'all_plasmids.fasta':
            continue
        seq_id, sites = scan_file(fpath)
        if seq_id is None:
            continue
        classification = classify_sites(sites)
        results.append({
            'name': fname.replace('.fasta', '').replace('.fa', ''),
            'seq_id': seq_id,
            'total': classification['total'],
            'cd': classification['cd'],
            'dc': classification['dc'],
            'dc_proportion': round(classification['dc_proportion'] * 100, 1),
        })
    return results


# ---------------------------------------------------------------------------
# Comparison table printer
# ---------------------------------------------------------------------------

def print_comparison_table(results, title="Comparison Results"):
    """Print a formatted comparison table.

    Parameters
    ----------
    results : list of dict
        Each dict must have: name, old_total, old_cd, old_dc, new_total, new_cd, new_dc
    """
    header = ("%-18s %-9s %-7s %-7s %-9s %-7s %-7s  %s" %
              ("Plasmid", "old_total", "old_CD", "old_DC",
               "new_total", "new_CD", "new_DC", "DC change"))
    sep = "-" * len(header)
    print("\n" + title)
    print(sep)
    print(header)
    print(sep)

    total_old_dc = 0
    total_old_all = 0
    total_new_dc = 0
    total_new_all = 0

    for r in results:
        name = r['name'][:17]
        err = r.get('error', '')
        if err:
            print("%-18s ERROR: %s" % (name, err[:60]))
            continue
        old_dc_pct = r.get('old_dc_pct', 0.0)
        new_dc_pct = r.get('new_dc_pct', 0.0)
        dc_change = "+%.1f%%" % (new_dc_pct - old_dc_pct) if new_dc_pct != old_dc_pct else "same"
        print("%-18s %-9d %-7d %-7d %-9d %-7d %-7d  (DC: %.1f%% -> %.1f%%, %s)" %
              (name, r['old_total'], r['old_cd'], r['old_dc'],
               r['new_total'], r['new_cd'], r['new_dc'],
               old_dc_pct, new_dc_pct, dc_change))
        total_old_dc += r['old_dc']
        total_old_all += r['old_total']
        total_new_dc += r['new_dc']
        total_new_all += r['new_total']

    old_dc_pct_total = round(total_old_dc / total_old_all * 100, 1) if total_old_all > 0 else 0.0
    new_dc_pct_total = round(total_new_dc / total_new_all * 100, 1) if total_new_all > 0 else 0.0
    print(sep)
    print("%-18s %-9d %-7d %-7d %-9d %-7d %-7d  DC: %.1f%% -> %.1f%%" %
          ("TOTAL", total_old_all, total_old_all - total_old_dc, total_old_dc,
           total_new_all, total_new_all - total_new_dc, total_new_dc,
           old_dc_pct_total, new_dc_pct_total))


def print_v2_table(results, title="V2 Scanner Results"):
    """Print v2-only results table."""
    header = "%-22s %-7s %-6s %-6s %-6s" % ("Plasmid", "total", "CD", "DC", "DC%")
    sep = "-" * 55
    print("\n" + title)
    print(sep)
    print(header)
    print(sep)
    total_all = 0
    total_cd = 0
    total_dc = 0
    for r in results:
        name = r['name'][:21]
        print("%-22s %-7d %-6d %-6d %-5.1f%%" %
              (name, r['total'], r['cd'], r['dc'], r['dc_proportion']))
        total_all += r['total']
        total_cd += r['cd']
        total_dc += r['dc']
    dc_pct = round(total_dc / total_all * 100, 1) if total_all > 0 else 0.0
    print(sep)
    print("%-22s %-7d %-6d %-6d %-5.1f%%" %
          ("TOTAL", total_all, total_cd, total_dc, dc_pct))


# ===================================================================
# Main: run full comparison when invoked directly
# ===================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(description='pdif v2 integration - compare scanners')
    parser.add_argument('--cache-dir', default=None,
                        help='Plasmid cache directory')
    parser.add_argument('--output-json', default=None,
                        help='Save results to JSON file')
    parser.add_argument('--mode', choices=['compare', 'v2-only', 'all'],
                        default='compare',
                        help='Comparison mode')
    parser.add_argument('--include-hall', action='store_true', default=True,
                        help='Include 4 Hall plasmids (download if needed)')
    args = parser.parse_args()

    cache_dir = args.cache_dir or _get_cache_dir()
    print("Cache dir: %s" % cache_dir)

    # Hall plasmid accessions
    hall_plasmids = {
        'p2ABAYE': 'NC_010402',
        'p2AB5075': 'CP008708',
        'pABV01': 'FM210331',
        'pAb242_9': 'KY984045',
        'pAb242_12': 'KY984046',
        'pAb242_25': 'KY984047',
    }

    comparison_results = []

    # 1. Download Hall plasmids
    print("\n--- Downloading Hall plasmids ---")
    fasta_paths = {}
    for name, accession in hall_plasmids.items():
        print("  %s (%s)..." % (name, accession))
        path = download_plasmid(accession, cache_dir)
        if path:
            fasta_paths[name] = path
            print("    -> %s" % path)
        else:
            print("    FAILED to download/cache")

    # 2. Run comparison
    if args.mode in ('compare', 'all'):
        print("\n--- Running comparison (old pdifFinder vs pwm scanner) ---")
        for name in ['pAb242_9', 'pAb242_12', 'pAb242_25',
                      'p2ABAYE', 'p2AB5075', 'pABV01']:
            if name not in fasta_paths:
                comparison_results.append({
                    'name': name, 'error': 'fasta not available',
                    'old_total': 0, 'old_cd': 0, 'old_dc': 0, 'old_dc_pct': 0.0,
                    'new_total': 0, 'new_cd': 0, 'new_dc': 0, 'new_dc_pct': 0.0,
                })
                continue
            print("  %s ..." % name)
            r = compare_plasmid(name, fasta_paths[name])
            comparison_results.append(r)
        print_comparison_table(comparison_results)

    # 3. V2-only analysis on ALL cached plasmids
    if args.mode in ('v2-only', 'all'):
        print("\n--- DC Proportion Analysis (pwm scanner, all cached plasmids) ---")
        all_results = analyze_all_cached_plasmids(cache_dir)
        print_v2_table(all_results)

        # Save DC proportion results
        if args.output_json:
            output = {
                'comparison': comparison_results,
                'v2_analysis': all_results,
                'total_plasmids_v2': len(all_results),
            }
            with open(args.output_json, 'w') as fh:
                json.dump(output, fh, indent=2)
            print("\nResults saved to %s" % args.output_json)

        # Print summary
        total_sites = sum(r['total'] for r in all_results)
        total_cd = sum(r['cd'] for r in all_results)
        total_dc = sum(r['dc'] for r in all_results)
        dc_pct = round(total_dc / total_sites * 100, 1) if total_sites > 0 else 0.0
        print("\nSUMMARY: %d plasmids, %d total sites, CD=%d, DC=%d (%.1f%%)" %
              (len(all_results), total_sites, total_cd, total_dc, dc_pct))


if __name__ == '__main__':
    main()
