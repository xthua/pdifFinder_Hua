#!/usr/bin/env python3
"""
Extract plasmid accessions from redundant.seed.fa, validate via NCBI eutils,
and produce tests/benchmark_data/plasmid_accessions.json.
"""

import json
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

# Allow unverified SSL context for NCBI eutils API in containerized env
ssl_ctx = ssl._create_unverified_context()

# ── Paths ──────────────────────────────────────────────────────────
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTA_PATH = os.path.join(PROJECT_ROOT, "PdifFinder", "data", "redundant.seed.fa")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "tests", "benchmark_data")
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "plasmid_accessions.json")

# ── Step 1: Parse FASTA ───────────────────────────────────────────
def parse_fasta(filepath):
    """Parse FASTA returning list of (header, sequence) tuples."""
    entries = []
    with open(filepath) as f:
        lines = [line.strip() for line in f if line.strip()]

    # Must be even number of lines (header + seq pairs)
    if len(lines) % 2 != 0:
        print(f"WARNING: odd number of non-empty lines ({len(lines)}), might be malformed")

    for i in range(0, len(lines), 2):
        header_line = lines[i]
        seq_line = lines[i + 1] if i + 1 < len(lines) else ""
        if not header_line.startswith(">"):
            print(f"WARNING: line {i+1} doesn't start with '>', skipping: {header_line[:40]}")
            continue
        header = header_line[1:]  # strip '>'
        entries.append((header, seq_line))

    return entries


# ── Step 2: Classify accessions ───────────────────────────────────
# NCBI accession pattern: 2 uppercase letters + digits + . + version
NCBI_ACC_RE = re.compile(r"^([A-Z]{2}\d+\.\d+)_(\d+)$")


def classify_header(header):
    """
    Return (accession, position, is_ncbi).
    - If header ends with _N (N = digits) and prefix looks like NCBI accession → NCBI
    - Otherwise → non-standard (valid=False)
    """
    m = NCBI_ACC_RE.match(header)
    if m:
        return m.group(1), int(m.group(2)), True
    # Non-standard identifier
    return header, None, False


# ── Step 3: Validate via NCBI eutils ──────────────────────────────
def validate_ncbi(accession, retries=3):
    """Check if an NCBI nucleotide accession exists via esearch."""
    params = urllib.parse.urlencode({"db": "nucleotide", "term": accession, "retmode": "json"})
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}"

    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PdifFinder/1.0"})
            with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as resp:
                data = json.loads(resp.read())
            count = int(data.get("esearchresult", {}).get("count", "0"))
            return count > 0
        except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError,
                OSError) as e:
            print(f"  Attempt {attempt + 1}/{retries} failed for {accession}: {e}")
            if attempt < retries - 1:
                time.sleep(2 ** attempt)  # exponential backoff
            else:
                return None  # unknown after all retries


# ── Main ──────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("Extracting plasmid accessions from redundant.seed.fa")
    print("=" * 60)

    # Parse
    entries = parse_fasta(FASTA_PATH)
    print(f"\nParsed {len(entries)} FASTA entries")

    # Group by accession
    from collections import OrderedDict
    accessions_map = OrderedDict()  # accession -> {"pdif_sequences": [...], "is_ncbi": bool}

    for header, seq in entries:
        acc, pos, is_ncbi = classify_header(header)
        if acc not in accessions_map:
            accessions_map[acc] = {"pdif_sequences": [], "is_ncbi": is_ncbi}
        accessions_map[acc]["pdif_sequences"].append(seq)

    print(f"\nUnique accessions: {len(accessions_map)}")
    for acc, data in accessions_map.items():
        print(f"  {acc:30s}  pdifs={len(data['pdif_sequences']):2d}  NCBI={data['is_ncbi']}")

    ncbi_accessions = [a for a, d in accessions_map.items() if d["is_ncbi"]]
    non_standard_accessions = [a for a, d in accessions_map.items() if not d["is_ncbi"]]

    print(f"\nNCBI accessions ({len(ncbi_accessions)}):")
    for a in ncbi_accessions:
        print(f"  {a}")
    print(f"\nNon-standard ({len(non_standard_accessions)}):")
    for a in non_standard_accessions:
        print(f"  {a}")

    # Validate NCBI accessions
    print(f"\n{'=' * 60}")
    print("Validating NCBI accessions via eutils API...")
    print(f"{'=' * 60}")

    results = []
    valid_ncbi_count = 0

    for acc, data in accessions_map.items():
        entry = {
            "id": acc,
            "pdif_count": len(data["pdif_sequences"]),
            "pdif_sequences": data["pdif_sequences"],
        }

        if data["is_ncbi"]:
            print(f"  Validating {acc}...", end=" ")
            sys.stdout.flush()
            valid = validate_ncbi(acc)
            if valid is True:
                print(f"✓ VALID")
                entry["valid"] = True
                valid_ncbi_count += 1
            elif valid is False:
                print(f"✗ NOT FOUND")
                entry["valid"] = False
            else:
                print(f"? UNKNOWN (API error)")
                entry["valid"] = None
            time.sleep(0.35)  # Rate limit: ~3 req/sec max for NCBI
        else:
            entry["valid"] = False
            print(f"  {acc:30s}  non-standard (skipped NCBI validation)")

        results.append(entry)

    non_standard_count = len(non_standard_accessions)

    # Build output
    output = {
        "accessions": results,
        "total": len(results),
        "valid_ncbi": valid_ncbi_count,
        "non_standard": non_standard_count,
        "_metadata": {
            "source": "redundant.seed.fa",
            "total_entries": len(entries),
            "description": "Plasmid accessions extracted from PdifFinder benchmark seed data",
        },
    }

    # Write output
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"Results written to: {OUTPUT_PATH}")
    print(f"  Total unique accessions: {len(results)}")
    print(f"  Valid NCBI accessions:   {valid_ncbi_count}")
    print(f"  Non-standard:            {non_standard_count}")
    print(f"{'=' * 60}")

    # Quick validation
    assert len(results) >= 12, f"Expected ≥12 accessions, got {len(results)}"
    assert valid_ncbi_count >= 10, f"Expected ≥10 valid NCBI, got {valid_ncbi_count}"

    # Verify every entry has required fields
    for r in results:
        assert "id" in r, f"Missing 'id' in entry"
        assert "valid" in r, f"Missing 'valid' in {r['id']}"
        assert "pdif_count" in r, f"Missing 'pdif_count' in {r['id']}"
        assert "pdif_sequences" in r, f"Missing 'pdif_sequences' in {r['id']}"
        assert len(r["pdif_sequences"]) == r["pdif_count"], \
            f"pdif_count mismatch for {r['id']}: {len(r['pdif_sequences'])} vs {r['pdif_count']}"
        for seq in r["pdif_sequences"]:
            assert len(seq) == 28, f"pdif seq not 28bp in {r['id']}: '{seq}' ({len(seq)}bp)"

    print("\n✓ All assertions passed")


if __name__ == "__main__":
    main()
