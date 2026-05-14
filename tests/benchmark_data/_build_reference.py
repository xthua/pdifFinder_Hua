#!/usr/bin/env python3
"""
Build literature_reference.csv from three sources:
1. Cameranesi et al. (2018) — 17 pdif sites across 3 plasmids
2. Blackwell & Hall (2017) — 8 pdif sites from 1 plasmid
3. redundant.seed.fa accessions — find pdif seed sequences in cached plasmid FASTAs
"""

import csv
import json
import sys
from pathlib import Path
from Bio import SeqIO

BASE = Path(__file__).resolve().parent

# Corrected PMIDs (verified against PubMed)
CAMERANESI_PMID = 29434581   # Not 29434599 (plan error) and not 29441038 (raw CSV error)
BLACKWELL_PMID = 28533235    # Verified
SHAO_PMID = 36426893         # Shao et al. 2023 Briefings in Bioinformatics

BLACKWELL_DESC = (
    "Blackwell GA, Hall RM. The tet39 Determinant and the msrE-mphE Genes "
    "in Acinetobacter Plasmids Are Each Part of Discrete Modules Flanked by "
    "Inversely Oriented pdif (XerC-XerD) Sites. Antimicrob Agents Chemother. "
    "2017;61(8):e00780-17. DOI: 10.1128/AAC.00780-17. "
    "8 pdif sites extracted from GenBank annotation of KY617771.1 (pS30-1) "
    "as misc_recomb features with /note=\"dif site\"."
)

SHAO_DESC = (
    "Shao C, Chen Y, Zhao Y, Wang Y, Song H, He M, Chen W, Xie Z, Qin S, "
    "Hua X. PdifFinder: a tool for identifying plasmid-borne dif sites "
    "and their cognate recombinases. Brief Bioinform. 2023;24(1):bbac521. "
    "DOI: 10.1093/bib/bbac521."
)


REQUIRED_COLUMNS = [
    "plasmid_accession", "pdif_start", "pdif_end", "pdif_sequence",
    "strand", "orientation", "source_pmid", "source_description",
    "curator_notes", "curation_date"
]

TARGET_CSV = BASE / "literature_reference.csv"


def load_plasmid_sequence(accession: str) -> str:
    """Load plasmid sequence from cached FASTA. Tries both . and _ variants."""
    for variant in [accession, accession.replace('.', '_')]:
        fasta_path = BASE / "plasmids" / f"{variant}.fasta"
        if fasta_path.exists():
            record = SeqIO.read(str(fasta_path), "fasta")
            return str(record.seq).upper()
    raise FileNotFoundError(f"No cached FASTA for {accession}")


def verify_coordinate(plasmid_seq: str, start: int, end: int, expected: str) -> bool:
    """Verify that plasmid_seq[start-1:end] matches expected (case-insensitive)."""
    extracted = plasmid_seq[start-1:end].upper()
    return extracted == expected.upper()


def extract_cameranesi_entries() -> list[dict]:
    """Parse cameranesi2018.csv and convert to target schema."""
    raw_path = BASE / "literature_raw" / "cameranesi2018.csv"
    entries = []
    
    with open(raw_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            accession = row["plasmid_accession"].strip()
            start = int(row["pdif_start"])
            end = int(row["pdif_end"])
            seq = row["pdif_sequence"].strip()
            strand = row.get("strand", "").strip() or "+"
            orientation = row.get("orientation", "").strip()
            
            # Verify and correct: Cameranesi CSV stores pdif_sequence in canonical
            # orientation. When strand='-', the actual plasmid sequence at [start:end]
            # is the reverse complement. Use the forward-strand sequence.
            try:
                plasmid_seq = load_plasmid_sequence(accession)
                actual_seq = plasmid_seq[start-1:end].upper()
                if actual_seq == seq.upper():
                    pass
                elif actual_seq == _reverse_complement(seq.upper()):
                    seq = actual_seq
                else:
                    print(f"WARNING: Cameranesi verification FAILED for {accession} "
                          f"{start}-{end}: expected={seq}, found={actual_seq}")
                    continue
            except FileNotFoundError:
                print(f"WARNING: No cached plasmid for {accession}, skipping")
                continue
            
            entries.append({
                "plasmid_accession": accession,
                "pdif_start": start,
                "pdif_end": end,
                "pdif_sequence": seq,
                "strand": strand,
                "orientation": orientation,
                "source_pmid": CAMERANESI_PMID,
                "source_description": row["source_description"].strip(),
                "curator_notes": row.get("curator_notes", "").strip(),
                "curation_date": "2026-05-10",
            })
    
    return entries


def extract_blackwell_entries() -> list[dict]:
    """Parse blackwell2017.csv and convert to target schema."""
    raw_path = BASE / "literature_raw" / "blackwell2017.csv"
    entries = []
    
    with open(raw_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            accession = row["plasmid_accession"].strip()
            start = int(row["pdif_start"])
            end = int(row["pdif_end"])
            seq = row["pdif_sequence"].strip()
            strand = row.get("strand", "").strip() or "+"
            
            entry = {
                "plasmid_accession": accession,
                "pdif_start": start,
                "pdif_end": end,
                "pdif_sequence": seq,
                "strand": strand,
                "orientation": "",
                "source_pmid": BLACKWELL_PMID,
                "source_description": BLACKWELL_DESC,
                "curator_notes": "",
                "curation_date": "2026-05-10",
            }
            
            # Verify against plasmid sequence
            try:
                plasmid_seq = load_plasmid_sequence(accession)
                verified = verify_coordinate(plasmid_seq, start, end, seq)
                if not verified:
                    extracted = plasmid_seq[start-1:end].upper()
                    print(f"WARNING: Blackwell verification FAILED for {accession} "
                          f"{start}-{end}: expected={seq}, found={extracted}")
                    continue
            except FileNotFoundError:
                print(f"WARNING: No cached plasmid for {accession}, skipping")
                continue
            
            entries.append(entry)
    
    return entries


def find_seed_in_plasmid(plasmid_seq: str, seed_seq: str) -> list[tuple[int, int]]:
    """Find all exact-match positions of seed_seq in plasmid_seq. Returns (start, end) 1-based."""
    positions = []
    pos = plasmid_seq.find(seed_seq.upper())
    while pos != -1:
        start_1based = pos + 1
        end_1based = pos + len(seed_seq)
        positions.append((start_1based, end_1based))
        pos = plasmid_seq.find(seed_seq.upper(), pos + 1)
    return positions


def extract_seed_entries() -> list[dict]:
    """For accessions NOT in Cameranesi/Blackwell, find pdif seed sequences in plasmids."""
    accessions_path = BASE / "plasmid_accessions.json"
    with open(accessions_path) as fh:
        data = json.load(fh)
    
    # Accessions already covered by Cameranesi or Blackwell
    covered = {"KY984045.1", "KY984046.1", "KY984047.1", "KY617771.1"}
    
    entries = []
    
    for acc in data["accessions"]:
        acc_id = acc["id"]
        
        # Skip non-standard and already-covered accessions
        if not acc["valid"] or acc_id in covered:
            continue
        
        # Load plasmid sequence
        try:
            plasmid_seq = load_plasmid_sequence(acc_id)
        except FileNotFoundError:
            print(f"WARNING: No cached plasmid for {acc_id}, skipping")
            continue
        
        # For each pdif seed sequence, find it in the plasmid
        for seed_seq in acc["pdif_sequences"]:
            positions = find_seed_in_plasmid(plasmid_seq, seed_seq)
            
            if not positions:
                print(f"WARNING: Seed sequence not found in {acc_id}: {seed_seq}")
                continue
            
            for start, end in positions:
                # Verify the extracted sequence
                extracted = plasmid_seq[start-1:end].upper()
                if extracted != seed_seq.upper():
                    print(f"WARNING: Verification mismatch for {acc_id} "
                          f"{start}-{end}: {extracted} != {seed_seq}")
                    continue
                
                entry = {
                    "plasmid_accession": acc_id,
                    "pdif_start": start,
                    "pdif_end": end,
                    "pdif_sequence": seed_seq,
                    "strand": "+",
                    "orientation": "",
                    "source_pmid": SHAO_PMID,
                    "source_description": SHAO_DESC,
                    "curator_notes": (
                        f"Located via exact string match of pdif seed sequence "
                        f"(from redundant.seed.fa) in cached plasmid sequence. "
                        f"This pdif site is present in the pdifFinder seed database "
                        f"as accession {acc_id}."
                    ),
                    "curation_date": "2026-05-10",
                }
                entries.append(entry)
    
    return entries


def _reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def deduplicate_entries(entries: list[dict]) -> list[dict]:
    """Remove duplicates by (plasmid_accession, pdif_start, pdif_end) key."""
    seen = set()
    unique = []
    for e in entries:
        key = (e["plasmid_accession"], e["pdif_start"], e["pdif_end"])
        if key not in seen:
            seen.add(key)
            unique.append(e)
    return unique


def write_csv(entries: list[dict]):
    """Write entries to literature_reference.csv with header comments."""
    # Sort: by plasmid_accession, then pdif_start
    entries.sort(key=lambda e: (e["plasmid_accession"], e["pdif_start"]))
    
    comment_block = """# literature_reference.csv — Curated pdif site references from published literature
#
# Schema:
#   plasmid_accession : str — NCBI GenBank accession of the plasmid
#   pdif_start        : int — Start coordinate of the pdif site (1-indexed)
#   pdif_end          : int — End coordinate of the pdif site (inclusive, 1-indexed)
#   pdif_sequence     : str — 28 bp pdif sequence (XerC(11) + spacer(6) + XerD(11))
#   strand            : str — '+' or '-', DNA strand
#   orientation       : str — 'direct' or 'inverted', relative orientation
#   source_pmid       : int — PubMed ID of the source publication
#   source_description: str — Brief human-readable citation or note
#   curator_notes     : str — Free-text notes from the curator
#   curation_date     : str — Date of curation in YYYY-MM-DD format
#
# Constraints:
#   - pdif_start < pdif_end
#   - pdif_sequence must be exactly 28 bp, valid DNA (A/C/G/T only)
#   - Combined unique key: (plasmid_accession, pdif_start, pdif_end)
#
# Source PMIDs:
#   29434581 - Cameranesi et al. 2018 (Front Microbiol)
#   28533235 - Blackwell & Hall 2017 (Antimicrob Agents Chemother)
#   36426893 - Shao et al. 2023 pdifFinder paper (Brief Bioinform)
#""".rstrip()
    
    with open(TARGET_CSV, "w", newline="", encoding="utf-8") as fh:
        fh.write(comment_block + "\n")
        writer = csv.DictWriter(fh, fieldnames=REQUIRED_COLUMNS)
        writer.writeheader()
        writer.writerows(entries)
    
    print(f"Wrote {len(entries)} entries to {TARGET_CSV}")


def main():
    print("=== Building literature_reference.csv ===")
    
    # Step 1: Cameranesi
    print("\n[1/3] Extracting Cameranesi et al. (2018) entries...")
    cameranesi = extract_cameranesi_entries()
    print(f"  Got {len(cameranesi)} verified entries")
    
    # Step 2: Blackwell
    print("\n[2/3] Extracting Blackwell & Hall (2017) entries...")
    blackwell = extract_blackwell_entries()
    print(f"  Got {len(blackwell)} verified entries")
    
    # Step 3: Seed sequences
    print("\n[3/3] Finding seed pdif sequences in cached plasmids...")
    seed_entries = extract_seed_entries()
    print(f"  Got {len(seed_entries)} entries from seed matching")
    
    # Merge and deduplicate
    all_entries = cameranesi + blackwell + seed_entries
    all_entries = deduplicate_entries(all_entries)
    
    # Statistics
    unique_plasmids = set(e["plasmid_accession"] for e in all_entries)
    unique_pmids = set(e["source_pmid"] for e in all_entries)
    
    print(f"\n=== SUMMARY ===")
    print(f"Total entries: {len(all_entries)}")
    print(f"Unique plasmids: {len(unique_plasmids)} → {sorted(unique_plasmids)}")
    print(f"Unique PMIDs: {len(unique_pmids)} → {sorted(unique_pmids)}")
    
    # Per-source breakdown
    print(f"\nPer-source breakdown:")
    for pmid in sorted(unique_pmids):
        count = sum(1 for e in all_entries if e["source_pmid"] == pmid)
        label = {CAMERANESI_PMID: "Cameranesi 2018",
                 BLACKWELL_PMID: "Blackwell 2017",
                 SHAO_PMID: "Shao 2023 (pdifFinder)"}.get(pmid, f"PMID {pmid}")
        print(f"  {label}: {count} entries")
    
    # Write CSV
    write_csv(all_entries)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
