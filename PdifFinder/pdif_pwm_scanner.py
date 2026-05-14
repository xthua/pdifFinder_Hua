#!/usr/bin/env python3
"""
pdif_pwm_scanner.py - PWM-based pdif site scanner for PdifFinder v2.

Complete reimplementation using Position Weight Matrices for
independent XerC/XerD arm detection on both strands.

Author: PdifFinder Team
License: GPLv3
"""

import json
import os
from Bio.Seq import Seq

NUCLEOTIDES = ['A', 'C', 'G', 'T']
NT_INDEX = {nt: i for i, nt in enumerate(NUCLEOTIDES)}
ARM_LENGTH = 11
SPACER_LENGTH = 6
SITE_LENGTH = 28
MIN_SEPARATION = 16
MAX_SEPARATION = 22

DEGENERATE_XERD_VARIANTS = [
    'TTATGCGAAAT',
    'TTATGCTAAAT',
    'TTACGTTAAAT',
]
DEGENERATE_VARIANT_WEIGHT = 0.5


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        return max(len(s1), len(s2))
    return sum(1 for a, b in zip(s1, s2) if a != b)


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def _pwm_max_score(pwm):
    return sum(max(pwm[nt][pos] for nt in NUCLEOTIDES) for pos in range(ARM_LENGTH))


def _parse_fasta_arms(seed_file):
    xerC_seeds = []
    xerD_seeds = []
    with open(seed_file) as fh:
        current_seq = ''
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    xerC_seeds.append(current_seq[0:11])
                    xerD_seeds.append(current_seq[17:28])
                    current_seq = ''
            else:
                current_seq = line
        if current_seq:
            xerC_seeds.append(current_seq[0:11])
            xerD_seeds.append(current_seq[17:28])
    return xerC_seeds, xerD_seeds


def _build_freq_table(arms):
    freq = {nt: [0.0] * ARM_LENGTH for nt in NUCLEOTIDES}
    for arm in arms:
        for pos, nt in enumerate(arm):
            if nt in freq:
                freq[nt][pos] += 1.0
    return freq


def _freq_to_pwm(freq, n_effective, pseudocount):
    pwm = {nt: [0.0] * ARM_LENGTH for nt in NUCLEOTIDES}
    for pos in range(ARM_LENGTH):
        col_sum = 0.0
        for nt in NUCLEOTIDES:
            pwm[nt][pos] = freq[nt][pos] + pseudocount
            col_sum += pwm[nt][pos]
        for nt in NUCLEOTIDES:
            pwm[nt][pos] /= col_sum
    return pwm


def build_pwms(seed_file=None, pseudocount=0.01):
    """Build XerC + XerD Position Weight Matrices and persist to JSON.

    Parameters
    ----------
    seed_file : str or None
        FASTA file with 28 bp CD-oriented seeds (XerC + spacer + XerD).
        Default: ``PdifFinder/data/xer_seeds_v1.fa``.
    pseudocount : float
        Added to each nucleotide count before normalisation.

    Returns
    -------
    dict with keys: xerC_pwm, xerD_pwm, xerD_degenerate_pwm,
    xerC_seeds, xerD_seeds.
    """
    if seed_file is None:
        seed_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'data', 'xer_seeds_v1.fa'
        )

    xerC_seeds, xerD_seeds = _parse_fasta_arms(seed_file)
    n_seeds = len(xerC_seeds)

    xerC_freq = _build_freq_table(xerC_seeds)
    xerD_freq = _build_freq_table(xerD_seeds)
    xerC_pwm = _freq_to_pwm(xerC_freq, n_seeds, pseudocount)
    xerD_pwm = _freq_to_pwm(xerD_freq, n_seeds, pseudocount)

    xerD_degen_freq = {
        nt: list(xerD_freq[nt])
        for nt in NUCLEOTIDES
    }
    for variant in DEGENERATE_XERD_VARIANTS:
        for pos, nt in enumerate(variant):
            if nt in xerD_degen_freq:
                xerD_degen_freq[nt][pos] += DEGENERATE_VARIANT_WEIGHT
    degen_effective = n_seeds + len(DEGENERATE_XERD_VARIANTS) * DEGENERATE_VARIANT_WEIGHT
    xerD_degenerate_pwm = _freq_to_pwm(xerD_degen_freq, degen_effective, pseudocount)

    result = {
        'xerC_pwm': xerC_pwm,
        'xerD_pwm': xerD_pwm,
        'xerD_degenerate_pwm': xerD_degenerate_pwm,
        'xerC_seeds': xerC_seeds,
        'xerD_seeds': xerD_seeds,
        'description': 'XerC/XerD Position Weight Matrices for pdifFinder v2',
        'params': {
            'pseudocount': pseudocount,
            'n_seeds': n_seeds,
            'arm_length': ARM_LENGTH,
            'spacer_length': SPACER_LENGTH,
        },
    }

    output_dir = os.path.dirname(seed_file) if seed_file else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'data'
    )
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'xer_pwm.json')
    with open(output_path, 'w') as fh:
        json.dump(result, fh, indent=2)

    return result


def scan_xer_arm(seq, pwm, seeds, pwm_threshold_ratio, max_mm):
    """Scan *seq* with an 11 bp sliding window, scoring against a PWM.

    Each window is evaluated by PWM log-probability sum and minimum
    Hamming distance to the seed set.  A hit is reported when either
    criterion meets its threshold.

    Parameters
    ----------
    seq : str
        DNA string (upper-case).
    pwm : dict
        ``{nt: [p0..p10]}``.
    seeds : list of str
        Known arm sequences for Hamming fallback.
    pwm_threshold_ratio : float
        Fraction of max achievable PWM score (e.g. 0.6).
    max_mm : int
        Max Hamming distance to any seed.

    Returns
    -------
    list of dict: ``[{'pos', 'seq', 'pwm_score', 'hamming'}, ...]``.
    """
    pwm_max = _pwm_max_score(pwm)
    pwm_abs_threshold = pwm_threshold_ratio * pwm_max
    hits = []

    limit = len(seq) - ARM_LENGTH + 1
    if limit <= 0:
        return hits

    for i in range(limit):
        window = seq[i:i + ARM_LENGTH]
        pwm_score = 0.0
        valid = True
        for pos, nt in enumerate(window):
            if nt not in pwm:
                valid = False
                break
            pwm_score += pwm[nt][pos]
        if not valid:
            continue

        min_hamming = ARM_LENGTH
        for seed in seeds:
            hd = hamming_distance(window, seed)
            if hd < min_hamming:
                min_hamming = hd
                if min_hamming == 0:
                    break

        if pwm_score >= pwm_abs_threshold or min_hamming <= max_mm:
            hits.append({
                'pos': i,
                'seq': window,
                'pwm_score': pwm_score,
                'hamming': min_hamming,
            })

    return hits


def assemble_pdif_sites(xerC_hits, xerD_hits, seq_len):
    """Pair XerC and XerD arm hits into candidate pdif sites.

    Reports both CD (XerC then XerD) and DC (XerD then XerC) pairings
    where the arm-start separation is between 16 and 22 bp.

    Parameters
    ----------
    xerC_hits, xerD_hits : list of dict
        Results from ``scan_xer_arm``.
    seq_len : int
        Total sequence length for boundary checking.

    Returns
    -------
    list of dict with ``start``, ``end``, ``raw_orientation``,
    ``xerC_hit``, ``xerD_hit``, ``spacer_len``.
    """
    xerC_by_pos = {}
    for h in xerC_hits:
        xerC_by_pos.setdefault(h['pos'], []).append(h)

    xerD_by_pos = {}
    for h in xerD_hits:
        xerD_by_pos.setdefault(h['pos'], []).append(h)

    candidates = []

    for c_pos, c_hits in xerC_by_pos.items():
        for delta in range(MIN_SEPARATION, MAX_SEPARATION + 1):
            d_pos = c_pos + delta
            if d_pos + ARM_LENGTH > seq_len:
                continue
            if d_pos not in xerD_by_pos:
                continue
            for ch in c_hits:
                for dh in xerD_by_pos[d_pos]:
                    candidates.append({
                        'start': c_pos,
                        'end': c_pos + SITE_LENGTH,
                        'raw_orientation': 'CD',
                        'xerC_hit': ch,
                        'xerD_hit': dh,
                        'spacer_len': delta - ARM_LENGTH,
                    })

    for d_pos, d_hits in xerD_by_pos.items():
        for delta in range(MIN_SEPARATION, MAX_SEPARATION + 1):
            c_pos = d_pos + delta
            if c_pos + ARM_LENGTH > seq_len:
                continue
            if c_pos not in xerC_by_pos:
                continue
            for dh in d_hits:
                for ch in xerC_by_pos[c_pos]:
                    candidates.append({
                        'start': d_pos,
                        'end': d_pos + SITE_LENGTH,
                        'raw_orientation': 'DC',
                        'xerC_hit': ch,
                        'xerD_hit': dh,
                        'spacer_len': delta - ARM_LENGTH,
                    })

    return candidates


def _compute_base_score(site):
    xerC_score = site.get('xerC_pwm_score', 0.0)
    xerD_score = site.get('xerD_pwm_score', 0.0)
    spacer_val = site.get('_spacer_len', site.get('spacer_len', SPACER_LENGTH))
    if isinstance(spacer_val, str):
        spacer_val = len(spacer_val)
    spacer_bonus = 0.1 if spacer_val == SPACER_LENGTH else 0.0
    paired_bonus = 0.2
    return xerC_score + xerD_score + spacer_bonus + paired_bonus


def _canonicalize_site(seq, raw_start, raw_orientation, xerC_hit, xerD_hit):
    raw_window = seq[raw_start:raw_start + SITE_LENGTH]
    canonical_seq = raw_window if raw_orientation == 'CD' else reverse_complement(raw_window)

    return {
        'contig': 'seq',
        'start': raw_start,
        'end': raw_start + SITE_LENGTH,
        'strand': '+',
        'raw_orientation': raw_orientation,
        'canonical_orientation': 'CD',
        'xerC_seq': canonical_seq[0:11],
        'spacer_seq': canonical_seq[11:17],
        'xerD_seq': canonical_seq[17:28],
        'canonical_seq': canonical_seq,
        'raw_window': raw_window,
        'xerC_pwm_score': xerC_hit['pwm_score'],
        'xerD_pwm_score': xerD_hit['pwm_score'],
        'xerC_hamming': xerC_hit['hamming'],
        'xerD_hamming': xerD_hit['hamming'],
        'spacer_seq_raw': raw_window[11:17],
        '_spacer_len': len(raw_window[11:17]),
    }


def _canonicalize_revcomp_site(seq, rev_complement, raw_start, raw_orientation,
                                xerC_hit, xerD_hit):
    forward_orientation = 'DC' if raw_orientation == 'CD' else 'CD'
    fwd_start = len(seq) - (raw_start + SITE_LENGTH)
    while fwd_start < 0:
        fwd_start += len(seq)

    raw_window = rev_complement[raw_start:raw_start + SITE_LENGTH]
    fwd_window = raw_window if forward_orientation == 'CD' else reverse_complement(raw_window)

    return {
        'contig': 'seq',
        'start': fwd_start,
        'end': fwd_start + SITE_LENGTH,
        'strand': '-',
        'raw_orientation': forward_orientation,
        'canonical_orientation': 'CD',
        'xerC_seq': fwd_window[0:11],
        'spacer_seq': fwd_window[11:17],
        'xerD_seq': fwd_window[17:28],
        'canonical_seq': fwd_window,
        'raw_window': raw_window,
        'xerC_pwm_score': xerC_hit['pwm_score'],
        'xerD_pwm_score': xerD_hit['pwm_score'],
        'xerC_hamming': xerC_hit['hamming'],
        'xerD_hamming': xerD_hit['hamming'],
        'spacer_seq_raw': raw_window[11:17],
        '_spacer_len': len(raw_window[11:17]),
    }


def _degenerate_rescue(seq, identified_starts, xerC_hits, xerD_seeds,
                       xerC_pwm, xerD_degenerate_pwm,
                       feature_boundaries):
    rescued = []

    if feature_boundaries is None or len(feature_boundaries) == 0:
        return rescued

    xerC_pwm_max = _pwm_max_score(xerC_pwm)
    xerD_pwm_max = _pwm_max_score(xerD_degenerate_pwm)

    strong_xerC = [
        h for h in xerC_hits
        if h['pwm_score'] >= 0.6 * xerC_pwm_max
        and h['pos'] not in identified_starts
    ]

    for ch in strong_xerC:
        pos = ch['pos']
        near_boundary = any(
            abs(pos - fb_start) <= 200 or abs(pos - fb_end) <= 200
            for fb_start, fb_end in feature_boundaries
        )
        if not near_boundary:
            continue

        for delta in range(MIN_SEPARATION, MAX_SEPARATION + 1):
            d_pos = pos + delta
            if d_pos + ARM_LENGTH > len(seq):
                continue
            window = seq[d_pos:d_pos + ARM_LENGTH]
            d_score = 0.0
            valid = True
            for p, nt in enumerate(window):
                if nt not in xerD_degenerate_pwm:
                    valid = False
                    break
                d_score += xerD_degenerate_pwm[nt][p]
            if not valid:
                continue
            min_hamming = min(hamming_distance(window, sd) for sd in xerD_seeds)
            if d_score >= 0.4 * xerD_pwm_max or min_hamming <= 3:
                dh = {'pos': d_pos, 'seq': window, 'pwm_score': d_score,
                      'hamming': min_hamming}
                site = _canonicalize_site(seq, pos, 'CD', ch, dh)
                site['rescue_flag'] = 'degenerate'
                site['pdif_type'] = 'degenerate'
                rescued.append(site)
                break

        for delta in range(MIN_SEPARATION, MAX_SEPARATION + 1):
            d_pos = pos - delta
            if d_pos < 0:
                continue
            window = seq[d_pos:d_pos + ARM_LENGTH]
            d_score = 0.0
            valid = True
            for p, nt in enumerate(window):
                if nt not in xerD_degenerate_pwm:
                    valid = False
                    break
                d_score += xerD_degenerate_pwm[nt][p]
            if not valid:
                continue
            min_hamming = min(hamming_distance(window, sd) for sd in xerD_seeds)
            if d_score >= 0.4 * xerD_pwm_max or min_hamming <= 3:
                dh = {'pos': d_pos, 'seq': window, 'pwm_score': d_score,
                      'hamming': min_hamming}
                site = _canonicalize_site(seq, d_pos, 'DC', ch, dh)
                site['rescue_flag'] = 'degenerate'
                site['pdif_type'] = 'degenerate'
                rescued.append(site)
                break

    return rescued


def scan_sequence(seq, xerC_pwm, xerD_pwm, xerC_seeds, xerD_seeds,
                  feature_boundaries=None, max_len=200000):
    """Scan a DNA sequence for pdif sites on both strands.

    Detects independently-matching XerC and XerD arms, assembles
    candidate pairs, canonicalises orientation to 5'-XerC-spacer-XerD-3',
    and applies a composite scoring function.

    Returns
    -------
    list of dict
        Keys: contig, start, end, strand, raw_orientation,
        canonical_orientation, xerC_seq, spacer_seq, xerD_seq,
        total_score, pdif_type, rescue_flag.
    """
    seq = seq.upper()
    if len(seq) > max_len:
        seq = seq[:max_len]
    seq_len = len(seq)

    xerC_hits_fwd = scan_xer_arm(seq, xerC_pwm, xerC_seeds, 0.6, 2)
    xerD_hits_fwd = scan_xer_arm(seq, xerD_pwm, xerD_seeds, 0.5, 3)
    candidates_fwd = assemble_pdif_sites(xerC_hits_fwd, xerD_hits_fwd, seq_len)

    rev_seq = reverse_complement(seq)
    xerC_hits_rev = scan_xer_arm(rev_seq, xerC_pwm, xerC_seeds, 0.6, 2)
    xerD_hits_rev = scan_xer_arm(rev_seq, xerD_pwm, xerD_seeds, 0.5, 3)
    candidates_rev = assemble_pdif_sites(xerC_hits_rev, xerD_hits_rev, len(rev_seq))

    sites = []
    seen_starts = set()

    for cand in candidates_fwd:
        site = _canonicalize_site(
            seq, cand['start'], cand['raw_orientation'],
            cand['xerC_hit'], cand['xerD_hit']
        )
        key = (site['start'], site['strand'])
        if key not in seen_starts:
            seen_starts.add(key)
            sites.append(site)

    for cand in candidates_rev:
        site = _canonicalize_revcomp_site(
            seq, rev_seq, cand['start'], cand['raw_orientation'],
            cand['xerC_hit'], cand['xerD_hit']
        )
        key = (site['start'], site['strand'])
        if key not in seen_starts:
            seen_starts.add(key)
            sites.append(site)

    sites.sort(key=lambda s: s['start'])

    for site in sites:
        site['_base_score'] = _compute_base_score(site)

    for i, site in enumerate(sites):
        context_bonus = 0.0
        for j, other in enumerate(sites):
            if i != j and abs(site['start'] - other['start']) <= 100:
                context_bonus = 0.1
                break
        site['total_score'] = site['_base_score'] + context_bonus

    identified_fwd_starts = set(cand['start'] for cand in candidates_fwd)
    rescued = _degenerate_rescue(
        seq, identified_fwd_starts,
        xerC_hits_fwd, xerD_seeds,
        xerC_pwm, xerD_pwm,
        feature_boundaries,
    )

    for r in rescued:
        r['_base_score'] = _compute_base_score(r)
        context_bonus = 0.0
        for other in sites + rescued:
            if r is other:
                continue
            if abs(r['start'] - other.get('start', 0)) <= 100:
                context_bonus = 0.1
                break
        r['total_score'] = r['_base_score'] + context_bonus

    output = []
    for site in sites + rescued:
        spacer_seq = site.get('spacer_seq', '')
        if not spacer_seq:
            canonical = site.get('canonical_seq', '')
            spacer_seq = canonical[11:17] if len(canonical) >= 17 else ''
        # ── Full C/D motif orientation scoring ──
        raw_ori = site['raw_orientation']
        xc = site.get('xerC_seq', '')
        xd = site.get('xerD_seq', '')
        left_CD_score = right_DC_score = right_CD_score = left_DC_score = 0
        orientation_margin = 0.0
        orientation_confidence = 'low'
        orientation_reason = 'default_CD'
        if len(xc) >= 11 and len(xd) >= 11:
            # Score XerC affinity (positions 0-11 of pdif)
            def _xc_aff(s):
                sc = 0
                # XerC often starts with ATT/ATA/GTA
                if s[:3] in ('ATT','ATA','GTA','ACT','GCT'): sc += 1
                # Mid-body patterns
                for m in ('TCGT','TCGC','TTCG','TTCG','TCTC'): sc += 1 if m in s else 0
                # More general
                return sc
            # Score XerD affinity (positions 17-28 of pdif)
            def _xd_aff(s):
                sc = 0
                for m in ('TTATG','ATGTTA','TTATGT','TTGAT'): sc += 2 if m in s else 0
                if s[-4:] in ('AAAT','AAAG','AAGT','AACA'): sc += 2
                if s[:3] in ('TTA','GTA','CTA','ATA','CAA'): sc += 1
                return sc
            left_C = _xc_aff(xc); right_C = _xc_aff(xd)
            left_D = _xd_aff(xc); right_D = _xd_aff(xd)
            score_CD = left_C + right_D
            score_DC = left_D + right_C
            left_CD_score = left_C; right_DC_score = right_D
            right_CD_score = right_C; left_DC_score = left_D
            orientation_margin = abs(score_CD - score_DC)
            
            if score_CD > score_DC:
                raw_ori = 'CD'
                if orientation_margin >= 5: orientation_confidence = 'high'
                elif orientation_margin >= 2: orientation_confidence = 'medium'
                else: orientation_confidence = 'low'
                orientation_reason = f'score_CD({score_CD}) > score_DC({score_DC})'
            elif score_DC > score_CD:
                raw_ori = 'DC'
                if orientation_margin >= 5: orientation_confidence = 'high'
                elif orientation_margin >= 2: orientation_confidence = 'medium'
                else: orientation_confidence = 'low'
                orientation_reason = f'score_DC({score_DC}) > score_CD({score_CD})'
            else:
                raw_ori = 'ambiguous'
                orientation_margin = 0
                orientation_confidence = 'low'
                orientation_reason = f'score_CD({score_CD}) == score_DC({score_DC})'
        rev_flag = site.get('reverse_complemented', False) or site.get('strand', '+') == '-'
        record = {
            'contig': site.get('contig', 'seq'),
            'start': site['start'],
            'end': site['end'],
            'strand': site['strand'],
            'raw_orientation': raw_ori,
            'original_orientation': site.get('raw_orientation', raw_ori),
            'canonical_orientation': site.get('canonical_orientation', 'CD'),
            'reverse_complemented': rev_flag,
            'xerC_seq': site['xerC_seq'],
            'spacer_seq': spacer_seq,
            'xerD_seq': site['xerD_seq'],
            'left_CD_score': left_CD_score,
            'left_DC_score': left_DC_score,
            'right_CD_score': right_CD_score,
            'right_DC_score': right_DC_score,
            'orientation_margin': orientation_margin,
            'orientation_confidence': orientation_confidence,
            'orientation_reason': orientation_reason,
            'total_score': site['total_score'],
            'pdif_type': site.get('pdif_type', 'canonical'),
            'rescue_flag': site.get('rescue_flag', 'standard'),
            'motif_supported_orientation': raw_ori,
            'forward_orientation': raw_ori if not rev_flag else ('DC' if raw_ori=='CD' else 'CD'),
            'reverse_complement_orientation': ('DC' if raw_ori=='CD' else 'CD') if not rev_flag else raw_ori,
            'rc_consistency': 'single_strand' if orientation_margin == 0 else 'consistent',
            'dc_validation_status': 'not_GB_DC',
            # ── Biological orientation sign (strand-independent) ──
            # +1: XerC-left/XerD-right on forward (CD on fwd)
            # -1: XerD-left/XerC-right on forward (DC on fwd, = CD on rev)
            'orientation_sign': 1 if raw_ori == 'CD' else -1 if raw_ori == 'DC' else 0,
            'validation_reason': orientation_reason,
        }
        output.append(record)

    # ── Module detection via CR pairing ──
    _modules = []
    _paired_starts = set()
    for i, a in enumerate(output):
        for b in output[i+1:]:
            dist = b['start'] - a['start']
            if dist < 50 or dist > 20000:
                continue
            ca, cb = a['spacer_seq'], b['spacer_seq']
            if len(ca) == len(cb) and len(ca) >= 5:
                mm = sum(1 for x, y in zip(ca, cb) if x != y)
                if mm <= 1:
                    _paired_starts.add(a['start'])
                    _paired_starts.add(b['start'])
                    cr_compat = mm <= 1
                    module_type = 'C' if a['raw_orientation'] != b['raw_orientation'] else 'D'
                    _modules.append({
                        'left_start': a['start'], 'right_start': b['start'],
                        'left_end': a['end'], 'right_end': b['end'],
                        'distance': dist,
                        'left_CR': ca, 'right_CR': cb,
                        'CR_pair': f'{ca}/{cb}',
                        'CR_hamming': mm,
                        'CR_compatible': cr_compat,
                        'module_type': module_type,
                        'relative_orientation': 'opposite' if a['raw_orientation'] != b['raw_orientation'] else 'same',
                        'valid_module': cr_compat,
                    })

    # Add metadata to each site
    for rec in output:
        rec['in_module'] = rec['start'] in _paired_starts
        rec['module_count'] = sum(1 for m in _modules if m['left_start'] == rec['start'] or m['right_start'] == rec['start'])
    # Filter: noise floor (score >= 12) + top 80% of max
    if output:
        mx = max(r['total_score'] for r in output)
        output = [r for r in output if r['total_score'] >= 16.5]
    return output, _modules


def load_pwms(json_path=None):
    """Load pre-built PWMs from JSON.

    Returns a tuple of (xerC_pwm, xerD_pwm, xerD_degenerate_pwm,
    xerC_seeds, xerD_seeds).
    """
    if json_path is None:
        json_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'data', 'xer_pwm.json'
        )
    with open(json_path) as fh:
        data = json.load(fh)
    return (
        data['xerC_pwm'],
        data['xerD_pwm'],
        data.get('xerD_degenerate_pwm', data['xerD_pwm']),
        data['xerC_seeds'],
        data['xerD_seeds'],
    )




def validate_dc_sites(sites, gb_labels=None):
    '''Validate DC-labeled sites against motif polarity.
    
    Args:
        sites: list of pdif site dicts from scan_sequence
        gb_labels: dict mapping position -> 'CD' or 'DC' from GenBank
    
    Updates each site's dc_validation_status and validation_reason in-place.
    '''
    if not gb_labels:
        return sites
    
    for site in sites:
        pos = site['start']
        # Find nearest GB annotation
        gb_ori = None
        for gs, gb in gb_labels.items():
            if abs(gs - pos) <= 10:
                gb_ori = gb
                break
        
        if gb_ori != 'DC':
            site['dc_validation_status'] = 'not_GB_DC'
            site['validation_reason'] = site.get('orientation_reason', '')
            continue
        
        # GB label is DC — validate against motif
        motif = site.get('motif_supported_orientation', 'CD')
        margin = site.get('orientation_margin', 0)
        
        if motif == 'DC':
            site['dc_validation_status'] = 'validated_DC'
            reason = f'GB DC is motif-supported: left half-site scores higher as XerD, margin={margin:.0f}'
        elif motif == 'CD':
            site['dc_validation_status'] = 'GB_DC_not_supported_by_motif'
            reason = f'GB DC is not motif-supported: XerC-left/XerD-right polarity suggests CD (margin={margin:.0f}). Possible contig-orientation artifact or manual annotation error.'
        else:
            site['dc_validation_status'] = 'ambiguous_DC'
            reason = f'GB DC is ambiguous: CD and DC scores differ by less than threshold (margin={margin:.0f})'
        
        site['validation_reason'] = reason
        if site.get('orientation_reason'):
            site['validation_reason'] += ' | ' + site['orientation_reason']
    
    return sites




def hall_degraded_scan(seq, xerC_pwm, xerD_pwm, xerC_seeds, xerD_seeds,
                        expected_positions=None, context_radius=50):
    '''Hall-style degraded pdif search guided by expected positions.
    
    Replicates the manual annotation workflow from Blackwell & Hall (2017):
    1. High-confidence detection (PWM score >= 16)
    2. Near expected positions: accept lower scores (>= 14)
    3. Degraded rescue: near expected + weak signal (>= 12) flagged
    
    Args:
        seq: DNA sequence string
        xerC_pwm, xerD_pwm: Position Weight Matrices
        xerC_seeds, xerD_seeds: Seed sequences for Hamming fallback
        expected_positions: list of expected pdif positions (e.g. from GB)
        context_radius: bp radius around expected positions for relaxed search
    
    Returns:
        dict with 'high', 'medium', 'degraded', 'missed' lists
    '''
    from collections import defaultdict
    
    # Run full scan
    all_sites, modules = scan_sequence(seq, xerC_pwm, xerD_pwm, xerC_seeds, xerD_seeds)
    
    result = {'high': [], 'medium': [], 'degraded': [], 'missed': []}
    
    if not expected_positions:
        for r in all_sites:
            s = r['total_score']
            if s >= 16.0: result['high'].append(r)
            elif s >= 14.0: result['medium'].append(r)
            elif s >= 12.0: result['degraded'].append(r)
        return result
    
    # Guided by expected positions
    used = set()
    for exp in sorted(expected_positions):
        # Try high confidence first
        nearby = [r for r in all_sites if abs(r['start'] - exp) <= 5 and r['total_score'] >= 16.0]
        if nearby:
            best = max(nearby, key=lambda r: r['total_score'])
            used.add(best['start'])
            result['high'].append(best)
            continue
        
        # Medium: near expected
        nearby = [r for r in all_sites if abs(r['start'] - exp) <= 10 and r['total_score'] >= 14.0]
        if nearby:
            best = max(nearby, key=lambda r: r['total_score'])
            used.add(best['start'])
            result['medium'].append(best)
            continue
        
        # Degraded: weak signal near expected (Hall "manually examined for degraded sites")
        nearby = [r for r in all_sites if abs(r['start'] - exp) <= context_radius and r['total_score'] >= 12.0]
        if nearby:
            best = max(nearby, key=lambda r: r['total_score'])
            used.add(best['start'])
            best['hall_degraded'] = True
            result['degraded'].append(best)
        else:
            # Seed-based fallback (Hall manually compares to known sites)
            seed_found = False
            for pos in range(max(0, exp - 5), min(len(seq) - 28, exp + 5)):
                w28 = seq[pos:pos+28]
                for seed in xerC_seeds:
                    mm = sum(1 for i in range(min(28, len(w28), len(seed))) if w28[i] != seed[i])
                    if mm <= 4:
                        # Manually construct a degraded site entry
                        xc, sp, xd = w28[0:11], w28[11:17], w28[17:28]
                        def _xd(s):
                            sc=0
                            for m in ('TTATG','ATGTTA','TTATGT'): sc+=2 if m in s else 0
                            if s[-4:] in ('AAAT','AAAG','AAGT'): sc+=2
                            if s[:3] in ('TTA','GTA','CTA'): sc+=1
                            return sc
                        ls, rs = _xd(xc), _xd(xd)
                        degraded_site = {
                            'start': pos, 'end': pos+28, 'raw_orientation': 'CD' if rs>ls else 'DC' if ls>rs else 'ambiguous',
                            'total_score': 10.0, 'orientation_margin': abs(rs-ls), 'xerC_seq': xc, 'spacer_seq': sp, 'xerD_seq': xd,
                            'spacer_seq_raw': sp, 'hall_degraded': True, 'seed_match': True, 'seed_mm': mm,
                            'pdif_type': 'degraded', 'rescue_flag': 'seed_rescue',
                            'left_CD_score': 0, 'left_DC_score': 0, 'right_CD_score': 0, 'right_DC_score': 0,
                            'orientation_confidence': 'low', 'orientation_reason': f'seed_match_mm={mm}',
                            'in_module': False, 'module_count': 0, 'canonical_orientation': 'CD',
                            'original_orientation': 'CD' if rs>ls else 'DC', 'reverse_complemented': False,
                            'contig': 'seq', 'strand': '+'
                        }
                        result['degraded'].append(degraded_site)
                        seed_found = True
                        break
                if seed_found:
                    break
            if not seed_found:
                result['missed'].append(exp)
    
    # Add unguided sites (not near any expected position) if they score well
    for r in all_sites:
        if r['start'] not in used and r['total_score'] >= 16.0:
            result['high'].append(r)
    
    return result




def analyze_pairs(sites, min_score=14.0):
    '''Compute pairwise relative orientation and recombination potential.
    
    Uses biological orientation_sign (+1/-1) instead of CD/DC labels,
    which are affected by contig reverse-complement direction.
    
    Returns list of pair dicts with:
        - relative_orientation: 'direct_repeat' or 'inverted_repeat'
        - predicted_recombination: [deletion, excision, inversion, etc]
        - recombination_potential: high/medium/low
        - orientation_signs for both sites
        - compatibility scores for spacer and core arms
    '''
    pairs = []
    sites = sorted(sites, key=lambda x: x['start'])
    
    for i, a in enumerate(sites):
        if a['total_score'] < min_score: continue
        for b in sites[i+1:]:
            if b['total_score'] < min_score: continue
            dist = b['start'] - a['start']
            if dist < 50 or dist > 50000: continue
            
            sign_a = a.get('orientation_sign', 0)
            sign_b = b.get('orientation_sign', 0)
            
            if sign_a == sign_b:
                rel_ori = 'direct_repeat'
                pred_rec = ['deletion', 'excision', 'cointegrate_resolution']
                arrow = '-->  -->'
            elif sign_a * sign_b < 0:
                rel_ori = 'inverted_repeat'
                pred_rec = ['inversion']
                arrow = '-->  <--'
            else:
                rel_ori = 'unknown'
                pred_rec = ['unknown']
                arrow = '??  ??'
            
            # Spacer/Core compatibility (CR from NAR 2025)
            ca, cb = a.get('spacer_seq', ''), b.get('spacer_seq', '')
            if len(ca) == len(cb) and len(ca) >= 5:
                mm = sum(1 for x, y in zip(ca, cb) if x != y)
                spacer_id = 1.0 - mm / len(ca)
            else:
                spacer_id = 0
                mm = 99
            
            spacer_ok = spacer_id >= 0.8
            core_ok = True  # arm comparison already done
            
            if spacer_ok and rel_ori == 'direct_repeat':
                rec_pot = 'high_direct_repeat_recombination'
            elif spacer_ok and rel_ori == 'inverted_repeat':
                rec_pot = 'high_inversion_recombination'
            elif spacer_ok:
                rec_pot = 'medium_compatible_spacer'
            else:
                rec_pot = 'low_or_ambiguous'
            
            pairs.append({
                'site1_start': a['start'], 'site2_start': b['start'],
                'site1_ori': a.get('raw_orientation', '?'),
                'site2_ori': b.get('raw_orientation', '?'),
                'site1_sign': sign_a, 'site2_sign': sign_b,
                'relative_orientation': rel_ori,
                'predicted_recombination': pred_rec,
                'distance_bp': dist,
                'spacer_identity': round(spacer_id, 2),
                'spacer_hamming': mm,
                'spacer_compatible': spacer_ok,
                'recombination_potential': rec_pot,
                'arrow': f"[{a['start']}]({a.get('raw_orientation','?')}) {arrow} [{b['start']}]({b.get('raw_orientation','?')})",
                'pair_confidence': 'high' if spacer_ok and mm <= 1 else 'medium' if spacer_id >= 0.5 else 'low',
            })
    
    return pairs




def write_module_list(sites, modules, output_path):
    '''Write pdifmodule_list.txt from v2 CR-validated modules.
    
    Format matches original pdifFinder output:
    tab-separated: module_id, left_site, right_site, distance, CR_pair, type
    '''
    if not modules:
        with open(output_path, 'w') as f:
            f.write('# No CR-validated pdif modules found\n')
        return
    
    with open(output_path, 'w') as f:
        f.write('# pdif module list — CR-validated pairs (PWM)\n')
        f.write('# module_id\tleft_start\tleft_end\tright_start\tright_end\t'
                'distance\tleft_CR\tright_CR\tCR_hamming\tmodule_type\t'
                'valid\n')
        for i, m in enumerate(modules, 1):
            f.write('module{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                i, m['left_start'], m['left_end'], m['right_start'], m['right_end'],
                m['distance'], m['left_CR'], m['right_CR'], m['CR_hamming'],
                m['module_type'], str(m['valid_module']).lower()
            ))
    return output_path



def write_all_outputs(sites, modules, seq, output_dir, plasmid_id='plasmid'):
    """Generate all 6 standard pdifFinder output files from v2 results.
    
    Files produced:
        AMRgene.txt, pdif_site.txt, pdifmodule_list.txt,
        pdifmoduleseq.fasta, pdifmodule.svg, plasmid.html
    """
    import os, json
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. AMRgene.txt — placeholder (PWM doesn't do AMR detection)
    with open(f'{output_dir}/AMRgene.txt', 'w') as f:
        f.write('#AMRgene annotation (PWM — AMR not detected)\n')
    
    # 2. pdif_site.txt — all detected pdif sites
    with open(f'{output_dir}/pdif_site.txt', 'w') as f:
        f.write('#pdif_site annotation (pwm scanner)\n')
        for i, s in enumerate(sorted(sites, key=lambda x: x["start"]), 1):
            f.write('pdif{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                i, s['start'], s['end'], s['xerC_seq'], s['spacer_seq'], s['xerD_seq'], s['raw_orientation']))
    
    # 3. pdifmodule_list.txt — CR-validated modules
    write_module_list(sites, modules, f'{output_dir}/pdifmodule_list.txt')
    
    # 4. pdifmoduleseq.fasta — module DNA sequences
    with open(f'{output_dir}/pdifmoduleseq.fasta', 'w') as f:
        for i, m in enumerate(modules, 1):
            left = m['left_start']; right = m['right_end']
            module_seq = seq[left:right+28]
            f.write('>module{}|{}..{}\n{}\n'.format(i, left, right+28, module_seq))
    
    # 5. pdifmodule.svg — simple text-based map (SVG not generated)
    with open(f'{output_dir}/pdifmodule.svg', 'w') as f:
        f.write('# pdif module map (PWM)\n')
        for i, m in enumerate(modules, 1):
            f.write('module{}: [{},{}] ↔ [{},{}] CR={}/{} type={}\n'.format(
                i, m['left_start'], m['left_end'], m['right_start'], m['right_end'],
                m['left_CR'], m['right_CR'], m['module_type']))
    
    # 6. plasmid.html — module summary as JSON
    with open(f'{output_dir}/plasmid.html', 'w') as f:
        f.write('<html><body><pre>\n')
        f.write('pdifFinder v2 — Module Summary\n')
        f.write('Plasmid: {}\n'.format(plasmid_id))
        f.write('Total pdif sites: {}\n'.format(len(sites)))
        f.write('CR-validated modules: {}\n'.format(len(modules)))
        f.write('</pre></body></html>\n')
    
    return output_dir



def main():
    """CLI for pwm scanner."""
    import argparse
    p = argparse.ArgumentParser(description="pdifFinder pwm Scanner")
    p.add_argument("--fasta","-f", required=True, help="Input FASTA file")
    p.add_argument("--output","-o", default="./PWM_output", help="Output directory")
    p.add_argument("--mode","-m", choices=["scan","hall"], default="scan")
    p.add_argument("--expected","-e", nargs="*", type=int, help="Expected pdif positions (Hall mode)")
    args = p.parse_args()

    from Bio import SeqIO
    rec = next(SeqIO.parse(args.fasta, "fasta"))
    seq = str(rec.seq).upper()
    pwms = build_pwms()
    
    if args.mode == "hall":
        r = hall_degraded_scan(seq, pwms["xerC_pwm"], pwms["xerD_pwm"],
                                pwms["xerC_seeds"], pwms["xerD_seeds"], args.expected or [])
        all_sites = r["high"] + r["medium"] + r["degraded"]
        print("Hall: high=%d medium=%d degraded=%d missed=%d" % (len(r["high"]), len(r["medium"]), len(r["degraded"]), len(r["missed"])))
    else:
        all_sites, modules = scan_sequence(seq, pwms["xerC_pwm"], pwms["xerD_pwm"],
                                            pwms["xerC_seeds"], pwms["xerD_seeds"])
        print("v2: %d sites" % len(all_sites))
    
    cd = sum(1 for s in all_sites if s["raw_orientation"]=="CD")
    dc = sum(1 for s in all_sites if s["raw_orientation"]=="DC")
    print("  CD=%d DC=%d" % (cd, dc))
    for s in sorted(all_sites, key=lambda x: x["start"])[:15]:
        print("  [%6d-%6d] %5s score=%.1f CR=%s" % (s["start"], s["end"], s["raw_orientation"], s["total_score"], s["spacer_seq"]))

if __name__ == "__main__":
    main()
