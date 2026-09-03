#!/usr/bin/env python3

import gzip
import random
import sys
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit(
        "Usage: normalizeQuery.py <query_fasta_or_fasta_gz> <random_seed>"
    )

input_path = Path(sys.argv[1])
random_seed = int(sys.argv[2])
canonical = set("ACGT")

iupac_replacements = {
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

random.seed(random_seed)
normalized_sequences = []
mapping = []
discarded_records = []
current_header = None
current_sequence = []

def get_output_filename(path: Path) -> str:
    filename = path.name.lower()

    suffixesWtGZ = [
        ".fasta.gz"
        ".fna.gz",
        ".fa.gz"
    ]

    for suffix in suffixesWtGZ:
        if filename.endswith(suffix):
            return path.name[:-len(suffix)] + ".q"

    return path.stem + ".q"

def flush_record():
    global current_header
    global current_sequence

    if current_header is None:
        return

    sequence = "".join(current_sequence).strip()

    if not sequence:
        discarded_records.append(
            (
                current_header,
                "Sequence is empty."
            )
        )

        return

    normalized = []
    changes = []
    invalid_characters = []

    for position, character in enumerate(sequence, start=1):
        upper = character.upper()

        if upper in canonical:
            replacement = upper

            if character != replacement:
                changes.append(
                    f"Position {position}: '{character}' replaced by "
                    f"'{replacement}'"
                )
        elif upper in set("RYSWKMBDHVN"):
            replacement = random.choice("ACGT")
            changes.append(
                f"Position {position}: IUPAC character '{character}' "
                f"replaced by '{replacement}'"
            )
        else:
            invalid_characters.append(
                f"Position {position}: invalid character '{character}'"
            )

        if not invalid_characters:
            normalized.append(replacement)

    if invalid_characters:
        discarded_records.append(
            (
                current_header,
                "Sequence discarded because it contains invalid characters: "
                + "; ".join(invalid_characters)
            )
        )
        return

    sequence_number = len(normalized_sequences) + 1
    normalized_sequence = "".join(normalized)

    normalized_sequences.append(normalized_sequence)
    mapping.append(
        (
            sequence_number,
            current_header,
            len(sequence),
            len(normalized_sequence),
            changes
        )
    )

open_function = gzip.open if input_path.name.endswith(".gz") else open

with open_function(input_path, "rt", encoding="utf-8") as handle:
    for raw_line in handle:
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith(">"):
            flush_record()
            current_header = line[1:].strip()
            current_sequence = []
        else:
            if current_header is None:
                raise SystemExit(
                    "The query file is not a valid FASTA file: "
                    "Found sequence data before first header."
                )

            current_sequence.append(line)

flush_record()

if not normalized_sequences:
    raise SystemExit(
        "No valid query sequences remain after query normalization."
    )

queryOutputFileName = get_output_filename(input_path)

with open(queryOutputFileName, "w", encoding="utf-8") as out:
    for sequence in normalized_sequences:
        out.write(sequence + '\n')

with open("query.map.tsv", "w", encoding="utf-8") as out:
    out.write("query_number\theader\toriginal_length\tnormalized_length\n")

    for number, header, original_length, normalized_length, _ in mapping:
        out.write(
            f"{number}\t{header}\t{original_length}\t{normalized_length}\n"
        )

with open("query.normalization.log", "w", encoding="utf-8") as out:
    out.write("PLAST query normalization\n")
    out.write(f"Input: {input_path}\n")
    out.write(f"Random seed: {random_seed}\n")
    out.write(f"Accepted queries: {len(normalized_sequences)}\n")
    out.write(f"Discarded queries: {len(discarded_records)}\n")

    changed_records = 0
    changed_characters = 0

    for number, header, _, _, changes in mapping:
        if changes:
            changed_records += 1
            changed_characters += len(changes)
            out.write(f"\nQuery {number}: {header}\n")

            for change in changes:
                out.write(change + "\n")

    for header, reason in discarded_records:
        out.write(f"\nDiscarded query: {header}\n")
        out.write(reason + "\n")

    if changed_records == 0:
        out.write("\nNo accepted sequence had to be changed.\n")
    else:
        out.write(f"\nChanged accepted sequences: {changed_records}\n")
        out.write(f"Total changes: {changed_characters}\n")
