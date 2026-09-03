#!/usr/bin/env python3

import re
import sys

if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: annotateResults.py <input.plast> <query.map.tsv> <output.plast>"
    )

query_headers = {}

with open(sys.argv[2], "r", encoding="utf-8") as mapping_handle:
    next(mapping_handle)

    for line in mapping_handle:
        line = line.rstrip("\n")

        if not line:
            continue

        fields = line.split("\t")

        if len(fields) < 2:
            continue

        query_number = fields[0]
        query_header = fields[1]

        query_headers[query_number] = query_header

query_pattern = re.compile(r"^Query\s+(\d+):\s*$")

with open(sys.argv[1], "r", encoding="utf-8") as input_handle, \
        open(sys.argv[3], "w", encoding="utf-8") as output_handle:

    for line in input_handle:
        match = query_pattern.match(line.rstrip("\n"))

        if match:
            query_number = match.group(1)
            query_header = query_headers.get(query_number)

            if query_header is not None:
                output_handle.write(f"Query {query_header}:\n")
                continue

        output_handle.write(line)
