#!/usr/bin/env python3
"""
build_seed_db.py — Build gold-standard Xer seed database from all available sources.

Pipeline:
  1. COLLECT — Gather all 28bp sequences from NAR2025, Cameranesi2018,
                 Blackwell2017, GenBank annotations, and existing seed DB.
  2. PARSE   — Split each into XerC(11) + CR(6) + XerD(11).
  3. ORIENT  — Normalise to CD orientation (XerC on left).
  4. DEDUP   — Remove duplicates by 28bp sequence.
  5. GROUP   — Group by CR sequence similarity.
  6. ALIGN   — Build per-group majority-rule consensus and PWM.
  7. BUILD   — Write xer_seeds_v1.fa, xer_seeds_v1.dc.fa,
                 xer_seeds_v1.pwm.json, build_report.json.
"""

import csv
import json
import os
import re
import sys
from collections import OrderedDict, defaultdict
from Bio.Seq import Seq

# ── Paths ───────────────────────────────────────────────────────────
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "PdifFinder", "data")
BENCHMARK_DIR = os.path.join(PROJECT_ROOT, "tests", "benchmark_data")
LITERATURE_RAW = os.path.join(BENCHMARK_DIR, "literature_raw")
PLASMIDS_DIR = os.path.join(BENCHMARK_DIR, "plasmids")
SCRIPTS_DIR = os.path.join(PROJECT_ROOT, "scripts")

# Output paths
OUT_CD = os.path.join(DATA_DIR, "xer_seeds_v1.fa")
OUT_DC = os.path.join(DATA_DIR, "xer_seeds_v1.dc.fa")
OUT_PWM = os.path.join(DATA_DIR, "xer_seeds_v1.pwm.json")
OUT_REPORT = os.path.join(BENCHMARK_DIR, "build_report.json")

# ── Validation status levels (higher = more trusted) ────────────────
STATUS_EXPERIMENTAL = 5   # NAR 2025 validated
STATUS_LITERATURE = 4     # Published literature (Cameranesi, Blackwell)
STATUS_CURATED = 3        # literature_reference.csv
STATUS_GENBANK = 2        # GenBank annotations
STATUS_SEEDDB = 1         # Existing seed database

# ── Orientation detection ───────────────────────────────────────────
# In CD orientation: XerC arm (positions 0:11) starts with AT...,
# XerD arm (positions 17:28) starts with TTAT... (highly conserved).
# We use TTAT at positions 17-20 as the primary CD-orientation signal.


def orient_to_cd(seq):
    """Ensure sequence is CD-oriented (XerC on left, XerD on right).

    Uses structural rules: XerD arms almost always start with 'TTAT'.
    If the right-side 11bp starts with 'TTAT', the sequence is CD.
    If the left-side 11bp starts with 'TTAT', flip to CD.
    Ambiguous cases are kept as-is.
    """
    left_start4 = seq[0:4]
    right_start4 = seq[17:21]

    left_is_xerd = (left_start4 == "TTAT")
    right_is_xerd = (right_start4 == "TTAT")

    # Clear CD signal: XerD on the right
    if right_is_xerd and not left_is_xerd:
        return seq, False

    # Clear DC signal: XerD on the left — flip
    if left_is_xerd and not right_is_xerd:
        return str(Seq(seq).reverse_complement()), True

    # Ambiguous: use first-base heuristic (XerC usually starts with A)
    if seq[0] == "A" and seq[17] != "A":
        return seq, False
    if seq[17] == "A" and seq[0] != "A":
        return str(Seq(seq).reverse_complement()), True

    return seq, False


def hamming(s1, s2):
    """Hamming distance between two equal-length strings."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))
    return sum(1 for a, b in zip(s1, s2) if a.upper() != b.upper())


# ── Helpers ─────────────────────────────────────────────────────────
def parse_fasta(filepath):
    """Parse FASTA file returning list of (header, sequence) tuples."""
    entries = []
    with open(filepath) as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    for i in range(0, len(lines), 2):
        hdr = lines[i]
        seq = lines[i + 1] if i + 1 < len(lines) else ""
        if hdr.startswith(">"):
            entries.append((hdr[1:], seq))
    return entries


def read_csv_delim(filepath):
    """Read CSV with auto-detected delimiter (comma or tab), skipping comment lines."""
    with open(filepath) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    if not lines:
        return []
    sample = "\n".join(lines[:20])
    dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
    reader = csv.DictReader(lines, dialect=dialect)
    rows = list(reader)
    return rows


# ── Step 1: COLLECT ─────────────────────────────────────────────────
def collect_nar2025():
    """Source A: NAR 2025 experimentally validated xrs sequences."""
    # XerC=11, CR=6, XerD=11 — already CD-oriented
    seqs = [
        ("xrs38_NAR2025", "ACTTCGTATAATCGCCATTATGTTAAAT",
         STATUS_EXPERIMENTAL, "NAR2025"),
        ("xrs105_NAR2025", "ATTTCGTATAAGGTGTATTATGTTAAGT",
         STATUS_EXPERIMENTAL, "NAR2025"),
        # xrs66 is identical to xrs38 — will be deduplicated
    ]
    return seqs


def collect_cameranesi2018():
    """Source B: Cameranesi 2018 pAb242 sites (17 sequences from CSV)."""
    csv_path = os.path.join(LITERATURE_RAW, "cameranesi2018.csv")
    rows = read_csv_delim(csv_path)
    seqs = []
    for row in rows:
        raw_seq = row.get("pdif_sequence", "").strip().upper()
        if len(raw_seq) != 28:
            continue
        name = row.get("pdif_site_name", "unknown").strip()
        seqs.append(
            ("%s_Cam2018" % name, raw_seq, STATUS_LITERATURE, "Cameranesi2018")
        )
    return seqs


def collect_blackwell2017():
    """Source C: Blackwell & Hall 2017 pS30-1 sites (8 sequences from CSV)."""
    csv_path = os.path.join(LITERATURE_RAW, "blackwell2017.csv")
    rows = read_csv_delim(csv_path)
    seqs = []
    for row in rows:
        raw_seq = row.get("pdif_sequence", "").strip().upper()
        if len(raw_seq) != 28:
            continue
        # Use accession + start as unique name
        acc = row.get("plasmid_accession", "unknown").strip()
        start = row.get("pdif_start", "0").strip()
        name = "%s_%s_BH2017" % (acc, start)
        seqs.append((name, raw_seq, STATUS_LITERATURE, "Blackwell2017"))
    return seqs


def collect_genbank():
    """Source D: GenBank annotated files — extract XerC/XerD features."""
    # Using BioPython SeqIO to parse GenBank features properly
    try:
        from Bio import SeqIO
    except ImportError:
        print("WARNING: Bio.SeqIO not available, skipping GenBank feature extraction")
        return []

    gb_files = [
        os.path.join(PLASMIDS_DIR, "KY617771_1.gb"),
        os.path.join(PLASMIDS_DIR, "KY984045_1.gb"),
        os.path.join(PLASMIDS_DIR, "KY984046_1.gb"),
    ]

    seqs = []
    for gb_path in gb_files:
        if not os.path.exists(gb_path):
            continue
        try:
            for record in SeqIO.parse(gb_path, "genbank"):
                accession = record.id
                for feat in record.features:
                    # Look for XerC/XerD in note= or label= qualifiers
                    notes = ""
                    if "note" in feat.qualifiers:
                        notes = " ".join(feat.qualifiers["note"])
                    if "label" in feat.qualifiers:
                        notes += " " + " ".join(feat.qualifiers["label"])
                    if "XerC" not in notes and "XerD" not in notes:
                        continue

                    # Extract the feature sequence
                    try:
                        feat_seq = str(feat.extract(record.seq)).upper()
                    except Exception:
                        continue

                    if len(feat_seq) != 28:
                        # Try to trim or skip
                        if len(feat_seq) < 20 or len(feat_seq) > 35:
                            continue
                        # Could be slightly longer — try to find 28bp window
                        # For now, skip non-28bp
                        continue

                    # Determine strand information
                    strand = "+"
                    if feat.location.strand == -1:
                        strand = "-"

                    loc_str = str(feat.location).replace(" ", "")
                    name = "%s_%s_GB" % (accession, loc_str)
                    seqs.append(
                        (name, feat_seq, STATUS_GENBANK, "GenBank")
                    )
        except Exception as e:
            print("WARNING: failed to parse %s: %s" % (gb_path, e))
            continue

    return seqs


def collect_seeddb():
    """Source E: Existing redundant.seed.fa (41 entries)."""
    fasta_path = os.path.join(DATA_DIR, "redundant.seed.fa")
    entries = parse_fasta(fasta_path)
    seqs = []
    for header, seq in entries:
        seq = seq.strip().upper()
        if len(seq) != 28:
            continue
        name = "%s_seedDB" % header
        seqs.append((name, seq, STATUS_SEEDDB, "SeedDB"))
    return seqs


def collect_literature_reference():
    """Additional source: curated literature_reference.csv."""
    ref_path = os.path.join(BENCHMARK_DIR, "literature_reference.csv")
    rows = read_csv_delim(ref_path)
    seqs = []
    for row in rows:
        raw_seq = row.get("pdif_sequence", "").strip().upper()
        if len(raw_seq) != 28:
            continue

        # Check if this is from the Cameranesi or Blackwell studies
        source_pmid = row.get("source_pmid", "").strip()
        src_tag = "LitRef"

        acc = row.get("plasmid_accession", "unknown").strip()
        start = row.get("pdif_start", "0").strip()
        strand = row.get("strand", "+").strip()

        name = "%s_%s_%s_%s" % (acc, start, strand, src_tag)
        seqs.append((name, raw_seq, STATUS_CURATED, "LiteratureReference"))
    return seqs


# ── Step 2: PARSE + ORIENT ──────────────────────────────────────────
def parse_28bp(seq):
    """Split 28bp into (xerc, cr, xerd). Returns None if invalid."""
    if len(seq) != 28:
        return None
    s = seq.upper()
    if not all(c in "ACGT" for c in s):
        return None
    return (s[0:11], s[11:17], s[17:28])


# ── Step 3: DEDUP ───────────────────────────────────────────────────
def deduplicate(entries):
    """
    Deduplicate by 28bp sequence, keeping the one with highest
    validation status. Ties broken by first-seen.

    entries: list of (name, seq_28bp, status, source_tag) tuples
    Returns: list of deduplicated (name, seq, status, source, flipped) tuples
    """
    best = OrderedDict()  # seq -> (name, status, source)
    raw_entries = list(entries)  # (name, seq, status, source)

    for name, seq, status, src in raw_entries:
        if seq not in best or status > best[seq][1]:
            best[seq] = (name, status, src)

    result = []
    for seq, (name, status, src) in best.items():
        result.append((name, seq, status, src, False))
    return result


# ── Step 4: GROUP by CR similarity ──────────────────────────────────
def group_by_cr(parsed_entries):
    """
    Group entries by CR sequence similarity.
    parsed_entries: list of (name, xerc, cr, xerd, seq, status, source, flipped)

    Returns: list of groups, each group is a list of parsed_entries.
    """
    # First pass: group by exact CR match
    cr_groups = defaultdict(list)
    for entry in parsed_entries:
        cr = entry[2]  # cr is at index 2
        cr_groups[cr].append(entry)

    # Second pass: merge groups where CR hamming distance <= 1
    groups = list(cr_groups.values())
    merged = True
    while merged:
        merged = False
        new_groups = []
        used = [False] * len(groups)
        for i in range(len(groups)):
            if used[i]:
                continue
            combined = list(groups[i])
            used[i] = True
            for j in range(i + 1, len(groups)):
                if used[j]:
                    continue
                close = False
                for e1 in combined:
                    cr1 = e1[2]
                    for e2 in groups[j]:
                        cr2 = e2[2]
                        if hamming(cr1, cr2) <= 1:
                            close = True
                            break
                    if close:
                        break
                if close:
                    combined.extend(groups[j])
                    used[j] = True
                    merged = True
            new_groups.append(combined)
        groups = new_groups

    return groups


# ── Step 5: ALIGN (build consensus) ─────────────────────────────────
def build_consensus(entries_in_group):
    """
    Build majority-rule consensus and conservation for a group.
    Returns (consensus_seq, conservation_list, pwm_dict).
    """
    seqs = [e[4] for e in entries_in_group]  # full 28bp
    if not seqs:
        return "", [], {}

    n = len(seqs)
    cols = [[s[i] for s in seqs] for i in range(28)]

    consensus = []
    conservation = []
    pwm = {"A": [], "C": [], "G": [], "T": []}

    for col in cols:
        counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        for base in col:
            if base in counts:
                counts[base] += 1

        # Majority rule
        max_base = max(counts, key=counts.get)
        max_count = counts[max_base]
        if max_count == 0:
            consensus.append("N")
            conservation.append(0.0)
        else:
            consensus.append(max_base)
            conservation.append(float(max_count) / n)

        # PWM frequencies (with pseudocount 0.1)
        pseudo = 0.1
        total = n + 4 * pseudo
        for base in "ACGT":
            pwm[base].append(round((counts[base] + pseudo) / total, 4))

    return "".join(consensus), conservation, pwm


# ── Step 6: BUILD outputs ───────────────────────────────────────────
def write_fasta(filepath, entries, dc=False):
    """Write FASTA file. entries: list of (name, seq, ...)"""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        for entry in entries:
            name = entry[0]
            seq = entry[4]  # full 28bp sequence
            if dc:
                seq = str(Seq(seq).reverse_complement())
            f.write(">%s\n" % name)
            f.write("%s\n" % seq)


def write_pwm(filepath, pwm_data):
    """Write PWM JSON file."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(pwm_data, f, indent=2)


def write_report(filepath, report_data):
    """Write build report JSON."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(report_data, f, indent=2)


# ── Main ─────────────────────────────────────────────────────────────
def main():
    print("=" * 70)
    print("  build_seed_db.py — Gold-standard Xer seed database builder")
    print("=" * 70)

    # ---- Step 1: COLLECT ----
    print("\n[Step 1] COLLECTING sequences from all sources...")

    sources = OrderedDict()
    sources["NAR2025"] = collect_nar2025()
    print("  NAR2025:              %d sequences" % len(sources["NAR2025"]))

    sources["Cameranesi2018"] = collect_cameranesi2018()
    print("  Cameranesi2018:       %d sequences" % len(sources["Cameranesi2018"]))

    sources["Blackwell2017"] = collect_blackwell2017()
    print("  Blackwell2017:        %d sequences" % len(sources["Blackwell2017"]))

    sources["GenBank"] = collect_genbank()
    print("  GenBank annotations:  %d sequences" % len(sources["GenBank"]))

    sources["SeedDB"] = collect_seeddb()
    print("  Existing seed DB:     %d sequences" % len(sources["SeedDB"]))

    sources["LiteratureRef"] = collect_literature_reference()
    print("  Literature reference: %d sequences" % len(sources["LiteratureRef"]))

    # Flatten all entries: (name, seq, status, source)
    all_entries = []
    for src_tag, entries in sources.items():
        all_entries.extend(entries)

    total_raw = len(all_entries)
    print("\n  Total raw entries: %d" % total_raw)

    # ---- Step 1b: DEDUP ----
    print("\n[Step 1b] DEDUPLICATING by 28bp sequence...")
    deduped = deduplicate(all_entries)
    print("  After dedup: %d unique sequences" % len(deduped))

    # ---- Step 2: PARSE + ORIENT ----
    print("\n[Step 2-3] PARSING and ORIENTING to CD...")
    parsed = []  # (name, xerc, cr, xerd, seq, status, source, flipped)
    skipped = 0
    flipped_count = 0
    for name, seq, status, src, _was_flipped in deduped:
        parts = parse_28bp(seq)
        if parts is None:
            skipped += 1
            continue
        xerc, cr, xerd = parts

        # Orient to CD
        oriented_seq, was_flipped = orient_to_cd(seq)
        if was_flipped:
            flipped_count += 1
            parts2 = parse_28bp(oriented_seq)
            if parts2 is None:
                skipped += 1
                continue
            xerc, cr, xerd = parts2

        parsed.append((name, xerc, cr, xerd, oriented_seq, status, src, was_flipped))

    print("  Parsed successfully: %d" % len(parsed))
    print("  Skipped (bad length/chars): %d" % skipped)
    print("  Flipped to CD: %d" % flipped_count)

    # Re-dedup after orientation (flipping might create duplicates)
    dedup_map = OrderedDict()
    for entry in parsed:
        seq = entry[4]
        status = entry[5]
        if seq not in dedup_map or status > dedup_map[seq][5]:
            dedup_map[seq] = entry
    parsed = list(dedup_map.values())
    print("  After re-dedup (post-orientation): %d unique" % len(parsed))

    # RC-dedup: remove reverse-complement duplicates
    rc_removed = 0
    seen_seqs = OrderedDict()
    for entry in parsed:
        seq = entry[4]
        status = entry[5]
        rc = str(Seq(seq).reverse_complement())
        if seq in seen_seqs or rc in seen_seqs:
            rc_removed += 1
            continue
        seen_seqs[seq] = entry
    parsed = list(seen_seqs.values())
    print("  After RC-dedup: %d unique (removed %d RC duplicates)" %
          (len(parsed), rc_removed))

    # ---- Step 3: GROUP ----
    print("\n[Step 4] GROUPING by CR similarity...")
    groups = group_by_cr(parsed)
    print("  Groups formed: %d" % len(groups))
    for i, grp in enumerate(groups):
        crs = set(e[2] for e in grp)
        print("    Group %d: %d members, CRs=%s" % (i + 1, len(grp), sorted(crs)))

    # ---- Step 4: ALIGN + CONSENSUS ----
    print("\n[Step 5] BUILDING consensus and PWM...")
    consensus_entries = []
    pwm_data = {}
    group_consensus_info = []

    for i, grp in enumerate(groups):
        consensus_seq, conservation, pwm = build_consensus(grp)

        # Representative CR (most common in group)
        cr_counter = defaultdict(int)
        for e in grp:
            cr_counter[e[2]] += 1
        rep_cr = max(cr_counter, key=cr_counter.get)

        group_id = "group_%02d_CR_%s" % (i + 1, rep_cr)
        consensus_name = "consensus_%s" % group_id

        info = {
            "group_id": group_id,
            "consensus_name": consensus_name,
            "rep_cr": rep_cr,
            "member_count": len(grp),
            "consensus_28bp": consensus_seq,
            "conservation": [round(c, 3) for c in conservation],
            "member_names": [e[0] for e in grp],
            "member_sequences": [e[4] for e in grp],
            "crs_in_group": sorted(set(e[2] for e in grp)),
        }
        group_consensus_info.append(info)

        # Add consensus to entries for FASTA output
        # Use the consensus as a new "seed" with CURATED status
        consensus_entries.append(
            (consensus_name, "", rep_cr, "",
             consensus_seq, STATUS_CURATED, "Consensus", False)
        )

        # Add PWM data
        pwm_data[group_id] = {
            "cr": rep_cr,
            "member_count": len(grp),
            "consensus": consensus_seq,
            "conservation": [round(c, 3) for c in conservation],
            "pwm": {base: [round(v, 4) for v in freqs] for base, freqs in pwm.items()},
        }

    # Combine unique seeds + consensus into final seed set
    all_seeds = list(parsed) + consensus_entries
    # One more dedup
    final_map = OrderedDict()
    for entry in all_seeds:
        seq = entry[4]
        if seq not in final_map:
            final_map[seq] = entry
    all_seeds = list(final_map.values())

    # ---- Step 5: BUILD output files ----
    print("\n[Step 6] WRITING output files...")

    # Canonicalize names for FASTA output
    fasta_entries = []
    seen_names = set()
    for entry in all_seeds:
        seq = entry[4]
        raw_name = entry[0]

        # Ensure unique names
        base_name = raw_name
        counter = 1
        while base_name in seen_names:
            base_name = "%s_v%d" % (raw_name, counter)
            counter += 1
        seen_names.add(base_name)

        fasta_entries.append((base_name,) + entry[1:])

    # Write CD FASTA
    write_fasta(OUT_CD, fasta_entries, dc=False)
    print("  CD seeds: %s (%d entries)" % (OUT_CD, len(fasta_entries)))

    # Write DC FASTA
    write_fasta(OUT_DC, fasta_entries, dc=True)
    print("  DC seeds: %s (%d entries)" % (OUT_DC, len(fasta_entries)))

    # Recalc entry list with DC names
    dc_fasta_entries = []
    for entry in fasta_entries:
        dc_fasta_entries.append(
            ("%s_DC" % entry[0],) + entry[1:]
        )

    # Write PWM JSON
    write_pwm(OUT_PWM, {
        "description": "Position Weight Matrix for Xer seed groups",
        "total_groups": len(groups),
        "total_seeds": len(parsed),
        "total_consensus": len(consensus_entries),
        "groups": pwm_data,
    })
    print("  PWM:      %s" % OUT_PWM)

    # ---- Build report ----
    src_counts = {}
    for entry in parsed:
        src = entry[6]
        src_counts[src] = src_counts.get(src, 0) + 1

    report = {
        "build_date": "2026-05-13",
        "pipeline": "build_seed_db.py",
        "summary": {
            "total_raw_sequences": total_raw,
            "after_dedup": len(deduped),
            "after_parse_orient": len(parsed),
            "total_final_seeds": len(all_seeds),
            "total_groups": len(groups),
            "total_consensus": len(consensus_entries),
            "flipped_to_CD": flipped_count,
        },
        "sources": {
            "NAR2025": {"raw": len(sources["NAR2025"])},
            "Cameranesi2018": {"raw": len(sources["Cameranesi2018"])},
            "Blackwell2017": {"raw": len(sources["Blackwell2017"])},
            "GenBank": {"raw": len(sources["GenBank"])},
            "SeedDB": {"raw": len(sources["SeedDB"])},
            "LiteratureRef": {"raw": len(sources["LiteratureRef"])},
        },
        "final_source_breakdown": src_counts,
        "groups": group_consensus_info,
        "output_files": {
            "cd_seeds": OUT_CD,
            "dc_seeds": OUT_DC,
            "pwm": OUT_PWM,
            "report": OUT_REPORT,
        },
    }

    write_report(OUT_REPORT, report)
    print("  Report:   %s" % OUT_REPORT)

    # ---- Verification ----
    print("\n[VERIFICATION]")
    assert os.path.exists(OUT_CD), "CD FASTA missing"
    assert os.path.exists(OUT_DC), "DC FASTA missing"
    assert os.path.exists(OUT_PWM), "PWM missing"
    assert os.path.exists(OUT_REPORT), "Report missing"

    # Check FASTA content
    with open(OUT_CD) as f:
        cd_lines = [l.strip() for l in f if l.strip()]
    cd_entries = len(cd_lines) // 2
    print("  CD FASTA entries: %d" % cd_entries)
    assert cd_entries == len(all_seeds), "CD count mismatch: %d vs %d" % (
        cd_entries, len(all_seeds)
    )

    with open(OUT_DC) as f:
        dc_lines = [l.strip() for l in f if l.strip()]
    dc_entries = len(dc_lines) // 2
    print("  DC FASTA entries: %d" % dc_entries)
    assert dc_entries == len(all_seeds), "DC count mismatch: %d vs %d" % (
        dc_entries, len(all_seeds)
    )

    # Verify all sequences are 28bp
    for i in range(0, len(cd_lines), 2):
        seq = cd_lines[i + 1]
        assert len(seq) == 28, "Entry %d: seq length %d != 28" % (i // 2 + 1, len(seq))
        assert all(c in "ACGT" for c in seq), "Entry %d: invalid chars in %s" % (
            i // 2 + 1, seq
        )

    # Check DC is reverse complement of CD
    for i in range(0, len(cd_lines), 2):
        cd_seq = cd_lines[i + 1]
        dc_seq = dc_lines[i + 1]
        expected_rc = str(Seq(cd_seq).reverse_complement())
        assert dc_seq == expected_rc, "DC mismatch at entry %d" % (i // 2 + 1)

    # Verify PWM JSON
    with open(OUT_PWM) as f:
        pwm_check = json.load(f)
    assert "groups" in pwm_check
    assert len(pwm_check["groups"]) == len(groups)

    print("\n" + "=" * 70)
    print("  BUILD COMPLETE — all outputs verified")
    print("  CD seeds:  %d" % cd_entries)
    print("  Groups:    %d" % len(groups))
    print("  Consensus: %d" % len(consensus_entries))
    print("=" * 70)


if __name__ == "__main__":
    main()
