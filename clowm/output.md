# Output

Typical workflow output includes the following files:

| File or Directory | Description |
|---|---|
| *plast_results/* | Directory containing final PLAST search results. |
| *plast_results/results.plast* | Original PLAST result file. |
| *plast_results/results.annotated.plast* | PLAST result file annotated with original FASTA query headers. |
| *plast_results/plast.stderr.txt* | Standard error output from PLAST Search. |
| *plast_results/search.command.txt* | PLAST Search command used by the workflow. |
| *result_stats.txt* | Summary statistics generated from the PLAST result file. |
| *graph/* | Graph files created during a graph construction step. |
| *indexed_graph/* | Graph and index files created during an index construction step. |
| *build.stdout.txt* | Standard output from PLAST graph construction. |
| *build.stderr.txt* | Standard error output from PLAST graph construction. |
| *build.command.txt* | PLAST Build command used by the workflow. |
| *build_index.stdout.txt* | Standard output from PLAST index construction. |
| *build_index.stderr.txt* | Standard error output from PLAST index construction. |

## Annotated Result File

PLAST internally identifies normalized queries by their order in the input file.

For example, an unannotated PLAST result (as reported in file *plast_results/results.plast*) may contain:

```text
Query 1:
```

The workflow replaces this internal number with the original FASTA header.

For example, if the original query file contains:

```text
>unknown_sequence_001
ACGTACGTACGT
```

the annotated result (as reported in file *plast_results/results.annotated.plast*) contains:

```text
Query unknown_sequence_001:
```

## Format of the (Annotated) PLAST Result File

The files *plast_results/results.plast* (raw output file) and *plast_results/results.annotated.plast* (annotated output file) are plain-text files. They contain the PLAST search results for all query sequences in the order in which they occur in the provided query file.

Result files are organized into query sections. Each section starts with a line of the following form:

```text
Query <id>:
```

In the raw output file *<id>* is simply a number. In the annotated output file, it is the original FASTA header of the corresponding query sequence, e.g.:

```text
Query unknown_sequence_001:
```

### Search Progress Messages

After the query header, PLAST may report progress messages such as:

```text
Searching for seeds
Extending seeds
Performing gapped extension
```

These lines describe the search and alignment-extension stages performed by PLAST. They are informational and are not alignment results themselves.

### Alignment Records

Each reported alignment starts with a summary line:

```text
Score: <score> Length: <length> E-value: <e-value>
```

For example:

```text
Score: 990  Length: 990 E-value: 8.67844e-128
```

Fields have the following meaning.

- *Score*: Alignment score calculated by PLAST. Higher scores generally indicate better alignments.
- *Length*: Length of the reported alignment, including alignment columns containing gaps where applicable.
- *E-value*: Expectation value assigned to the alignment. Smaller values indicate stronger statistical significance.

One query may have several alignment records. The records are separated by their summary lines and are reported increasingly ordered by expectation value.

### Alignment Layout

An alignment is displayed using three lines for each alignment block:

```text
Query:    <start> <query sequence> <end>
          <match representation>
Graph:    <start> <graph sequence> <end>
```

For example:

```text
Query:    886   TGTCAAAAAAGCCTCCGAGCTGTTCCC-CCAAATCCAACAATTGACAAGGGATGGGTCTG    944
                ||||| ||||||||||| || | |||| |||   ||||  || || |  || ||||||||
Graph:    0     TGTCAGAAAAGCCTCCG-GCCGGTCCCACCATCACCAAAGATCGATAGAGGTTGGGTCTG    0
```

The three lines contain

- the query sequence and its coordinates,
- a character-by-character match representation, and
- the corresponding sequence from the pangenome graph.

The coordinate values refer to positions in the respective query sequence. A hyphen (-) represents a gap introduced by the alignment.

The match line uses characters to show the relationship between aligned positions. In particular, *|* indicates matching bases and spaces indicate mismatches or alignment positions without an identical base.

Long alignments are split into several consecutive alignment blocks for readability. The coordinate values allow the blocks to be placed in the context of the complete reported alignment.

> **Important:** Sequence coordinates are always 0 in lines representing the graph sequence. This is because a graph sequence is not a linear sequence and a linear coordinate system does not apply to it. 

### Color Coverage Information

If the Parameter `--report_colors`` is enabled, PLAST additionally reports the graph colors supporting different parts of an alignment.

These lines have the following form:

```text
Color set ending at alignment position <position> <color_1> <color_2> ...
```

For example:

```text
Color set ending at alignment position 639 (position 639 in the query sequence) genome_01.fa
Color set ending at alignment position 641 (position 639 in the query sequence) genome_01.fa genome_02.fasta
```

A color usually corresponds to one input sequence file used to construct the pangenome graph. The line states that the listed colors support the alignment up to the indicated alignment position. To facilitate mappings of color coverage to the query sequence an additional position is given in brackets. This position disregards gaps in the alignments and, thus, may be different from the previously stated position in the alignment.

Depending on the alignment and the selected graph, the output may contain

- one color,
- several colors, or
- several consecutive color coverage records.

If `--report_colors` is disabled, these lines are not written.

### Example of a Complete Result Section

```text
Query unknown_sequence_001:
Searching for seeds
Extending seeds
Performing gapped extension
Score: 990  Length: 990 E-value: 8.67844e-128
Query:    1 ATGTGGGACTCGTCATACATGCAACAAGTGAGTGAGGGACTGATGACTGGAAAAGTTCCA    60
            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Graph:    0 ATGTGGGACTCGTCATACATGCAACAAGTGAGTGAGGGACTGATGACTGGAAAAGTTCCA    0

...

Score: 38   Length: 104 E-value: 0.000939606
Query:    886   TGTCAAAAAAGCCTCCGAGCTGTTCCC-CCAAATCCAACAATTGACAAGGGATGGGTCTG    944
                ||||| ||||||||||| || | |||| |||   ||||  || || |  || ||||||||
Graph:    0     TGTCAGAAAAGCCTCCG-GCCGGTCCCACCATCACCAAAGATCGATAGAGGTTGGGTCTG    0

```

The next query section starts with the next *Query* line:

```text
Query unknown_sequence_002:
```

All subsequent alignment records belong to this query until the next query header or the end of the file is reached.

## Result Statistics

After a PLAST search, the workflow runs the script *showPLASTresStats.py* on the PLAST result file.

The generated summary is written to *result_stats.txt*.

This file is included in the final workflow output and can also be downloaded from the selected CloWM output bucket.
