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
            #Testing
            # print("Test 2: Suffix of compressed file detected.")

            return path.name[:-len(suffix)] + ".q"

    #Testing
    # print("Tests 3 and 4: No known suffix of compressed data detected.")

    return path.stem + ".q"

def flush_record():
    global current_header
    global current_sequence

    if current_header is None:
        #Testing
        # print("Tests 2, 3, 4, and 5: Variable current_header is empty")

        return

    #Testing
    # print("Tests 2, 3, 4, and 5: Variable current_header is not empty")

    sequence = "".join(current_sequence).strip()

    if not sequence:
        #Testing
        # print("Tests 4 and 5: FASTA record has an empty sequence")

        discarded_records.append(
            (
                current_header,
                "Sequence is empty."
            )
        )

        return

    #Testing
    # print("Tests 2, 3, 4, and 5: FASTA record has a non-empty sequence")

    normalized = []
    changes = []
    invalid_characters = []

    for position, character in enumerate(sequence, start=1):
        upper = character.upper()

        if upper in canonical:
            #Testing
            # print("Tests 2, 3, 4, and 5: Non-arbitrary nucleotid detected")

            replacement = upper

            if character != replacement:
                #Testing
                # print("Tests 2 and 3: Lowercase non-arbitrary character detected")

                changes.append(
                    f"Position {position}: '{character}' replaced by "
                    f"'{replacement}'"
                )
        elif upper in set("RYSWKMBDHVN"):
            #Testing
            # print("Tests 2 and 3: Arbitrary nucleotid detected")

            replacement = random.choice("ACGT")
            changes.append(
                f"Position {position}: IUPAC character '{character}' "
                f"replaced by '{replacement}'"
            )
        else:
            #Testing
            # print("Test 4 and 5: Found non-nucleotid character in sequence")

            invalid_characters.append(
                f"Position {position}: invalid character '{character}'"
            )

            # replacement = random.choice("ACGT")
            # changes.append(
            #     f"Position {position}: ungültiges Zeichen '{character}' "
            #     f"durch '{replacement}' ersetzt"
            # )

        # if character != replacement and upper in canonical:
        #     changes.append(
        #         f"Position {position}: '{character}' replaced by "
        #         f"'{replacement}' "
        #     )

        if not invalid_characters:
            #Testing
            # print("Tests 2, 3, 4, and 5: Found sequence with no invalid characters")

            normalized.append(replacement)

    if invalid_characters:
        #Testing
        # print("Test 4 and 5: Discard sequence, because it contains invalid characters")

        discarded_records.append(
            (
                current_header,
                "Sequence discarded because it contains invalid characters: "
                + "; ".join(invalid_characters)
            )
        )
        return

    #Testing
    # print("Tests 2, 3, and 4: Not discarded sequence")

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

#Testing
# if input_path.name.endswith(".gz"):
#     print("Test 2: Received compressed file")
# else:
#     print("Test 1, 3, 4, and 5: Received uncompressed file")

open_function = gzip.open if input_path.name.endswith(".gz") else open

with open_function(input_path, "rt", encoding="utf-8") as handle:

# with open_function(input_path, "rt") as handle:

    for raw_line in handle:
        line = raw_line.strip()

        if not line:
            #Testing
            # print("Tests 2, 3, 4, and 5: Found empty line in input file")

            continue

        if line.startswith(">"):
            #Testing
            # print("Tests 2, 3, 4, and 5: Found beginning of FASTA header")

            flush_record()
            current_header = line[1:].strip()
            current_sequence = []
        else:
            #Testing
            # print("Tests 1, 2, 3, 4, and 5: Found a line that does not start with a FASTA header")

            if current_header is None:
                #Testing
                # print("Test 1: Found a line with content without to have seen a FASTA header before")

                raise SystemExit(
                    "The query file is not a valid FASTA file: "
                    "Found sequence data before first header."
                )

            current_sequence.append(line)

flush_record()

if not normalized_sequences:
    #Testing
    # print("Test 5: No valid sequences remain after query normalization")

    raise SystemExit(
        "No valid query sequences remain after query normalization."
    )

#Testing
# print("Tests 2, 3, and 4: Some valid sequences remain after query normalization")

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
            #Testing
            # print("Tests 2 and 3: Changes had to be made during normalization")

            changed_records += 1
            changed_characters += len(changes)
            out.write(f"\nQuery {number}: {header}\n")

            for change in changes:
                out.write(change + "\n")

    for header, reason in discarded_records:
        out.write(f"\nDiscarded query: {header}\n")
        out.write(reason + "\n")

    if changed_records == 0:
        #Testing
        # print("Test 4: None of the sequences that were accepted had to be changed")

        out.write("\nNo accepted sequence had to be changed.\n")
    else:
        #Testing
        # print("Tests 2 and 3: An accepted sequence had to be changed")

        out.write(f"\nChanged accepted sequences: {changed_records}\n")
        out.write(f"Total changes: {changed_characters}\n")
