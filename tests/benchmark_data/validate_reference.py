#!/usr/bin/env python3
"""
validate_reference.py — Validate literature_reference.csv against schema constraints.

Exits with code 0 if valid, 1 if any errors are found.

Uses only the Python standard library (no BioPython dependency).
"""

import csv
import re
import sys
from pathlib import Path

REQUIRED_COLUMNS = [
    "plasmid_accession",
    "pdif_start",
    "pdif_end",
    "pdif_sequence",
    "strand",
    "orientation",
    "source_pmid",
    "source_description",
    "curator_notes",
    "curation_date",
]

VALID_DNA_RE = re.compile(r"^[ACGTacgt]+$")
DEFAULT_CSV_PATH = Path(__file__).resolve().parent / "literature_reference.csv"


def _non_comment_lines(fh):
    """Yield non-comment lines (skip lines starting with '#')."""
    for line in fh:
        stripped = line.lstrip()
        if stripped and not stripped.startswith("#"):
            yield line


def validate(csv_path: str | Path) -> bool:
    """Validate the reference CSV. Returns True if valid, False otherwise."""
    errors: list[str] = []
    seen_keys: set[tuple[str, int, int]] = set()
    row_count = 0

    path = Path(csv_path)

    if not path.exists():
        print(f"ERROR: File not found: {path}")
        return False

    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(_non_comment_lines(fh))
        header = reader.fieldnames

        if header is None:
            print("ERROR: CSV file is empty or has no header row.")
            return False

        missing_cols = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing_cols:
            print(
                f"ERROR: Missing required columns: {', '.join(missing_cols)}"
            )
            return False

        extra_cols = [c for c in header if c not in REQUIRED_COLUMNS]
        if extra_cols:
            print(
                f"NOTE: Extra columns found (not in schema): {', '.join(extra_cols)}"
            )

        for row in reader:
            row_count += 1
            line_num = row_count + 1

            accession = row.get("plasmid_accession", "").strip()
            start_str = row.get("pdif_start", "").strip()
            end_str = row.get("pdif_end", "").strip()
            seq = row.get("pdif_sequence", "").strip()
            source_pmid_str = row.get("source_pmid", "").strip()

            if not accession:
                errors.append(
                    f"Line {line_num}: plasmid_accession is empty"
                )

            try:
                start = int(start_str)
            except (ValueError, TypeError):
                errors.append(
                    f"Line {line_num}: pdif_start '{start_str}' is not a valid integer"
                )
                start = None
            try:
                end = int(end_str)
            except (ValueError, TypeError):
                errors.append(
                    f"Line {line_num}: pdif_end '{end_str}' is not a valid integer"
                )
                end = None

            if start is not None and end is not None:
                if start >= end:
                    errors.append(
                        f"Line {line_num}: pdif_start ({start}) >= pdif_end ({end})"
                    )

                key = (accession, start, end)
                if key in seen_keys:
                    errors.append(
                        f"Line {line_num}: Duplicate entry — "
                        f"(accession={accession}, start={start}, end={end}) "
                        f"already seen"
                    )
                else:
                    seen_keys.add(key)

            if not seq:
                errors.append(
                    f"Line {line_num}: pdif_sequence is empty"
                )
            elif len(seq) != 28:
                errors.append(
                    f"Line {line_num}: pdif_sequence length is {len(seq)}, expected 28"
                )
            elif not VALID_DNA_RE.match(seq):
                invalid_chars = set(seq.upper()) - {"A", "C", "G", "T"}
                errors.append(
                    f"Line {line_num}: pdif_sequence contains invalid DNA "
                    f"characters: {''.join(sorted(invalid_chars))}"
                )

            if not source_pmid_str:
                errors.append(
                    f"Line {line_num}: source_pmid is empty"
                )
            else:
                try:
                    pmid = int(source_pmid_str)
                    if pmid <= 0:
                        errors.append(
                            f"Line {line_num}: source_pmid '{source_pmid_str}' "
                            f"must be a positive integer"
                        )
                except (ValueError, TypeError):
                    errors.append(
                        f"Line {line_num}: source_pmid '{source_pmid_str}' "
                        f"is not a valid integer"
                    )

    if row_count == 0:
        print("WARNING: CSV has a header but no data rows.")

    if errors:
        print(f"\nFound {len(errors)} validation error(s):")
        for err in errors:
            print(f"  {err}")
        return False
    else:
        print(f"All validations passed ({row_count} row(s) checked).")
        return True


def main() -> int:
    csv_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CSV_PATH
    valid = validate(csv_path)
    return 0 if valid else 1


if __name__ == "__main__":
    sys.exit(main())
