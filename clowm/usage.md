# Usage

PLAST searches for high-scoring local alignments between DNA query sequences and a pangenome represented as a compacted colored de Bruijn graph.

The workflow can either

- build a new pangenome graph from reference sequences and/or read
  sequences,
- build a PLAST index for an existing graph,
- search an existing graph, or
- build a graph and subsequently search it.

# Running PLAST on CloWM

## Uploading Input Data

Before starting the workflow, all input data must be uploaded to a [CloWM data bucket](https://cmg.pages.ub.uni-bielefeld.de/clowm/clowm-user-documentation/02_user_data_upload/).

To manage input data:

1. [Log in to CloWM](https://cmg.pages.ub.uni-bielefeld.de/clowm/clowm-user-documentation/00A_access_clowm).
2. Open **Files > My Data Buckets**.
3. Select an accessible bucket.
4. Optionally create folders to organize input data.
5. Upload files or directories using the **Upload File** button.

```text
initial-bucket-xyz/
├── graphs/
│   ├── myPangenome.gfa.gz
│   ├── myPangenome.color.bfg
│   └── myPangenome.idx
├── queries/
│   └── myQueries.fasta
├── references/
│   ├── genome_01.fasta
│   ├── genome_02.fasta
│   └── selected_genomes.fof
├── reads/
│   ├── sample_01.fastq.gz
│   └── sample_02.fastq.gz
└── search_sets/
    └── selected_colors.txt
```

## Starting a Workflow

To start PLAST on CloWM:

1. Click on the green **Launch vX.Y.Z** button in the top of this page to run the latest workflow version. To run a workflow version other than *vX.Y.Z*, choose a version from the drop-down menu below the green button.
2. Select the desired parameter view:
   - **Simple** for essential input parameters;
   - **Advanced** for additional settings;
   - **Expert** for all available parameters.
3. Select input files and directories from a CloWM data bucket.
4. Select an output bucket and output folder using Parameter `--outdir`.
5. Click **Launch**.

After starting the workflow, CloWM opens the workflow execution report. The
report shows the status of the workflow, individual workflow processes, and
resource usage.

# Workflow Modes

The Parameter `--mode` determines which workflow steps are performed.

Possible values are `Build`, `Search`, and `Build and Search`.

## Build

The Build mode creates a new pangenome graph or builds a PLAST index for an existing graph.

A query sequence file is not required in this mode.

### Build from Sequence Files

If sequence input directories or a ZIP archive are provided, the workflow builds

1. a compacted colored de Bruijn graph, and
2. a PLAST search index.

The workflow accepts

- reference sequences through Parameter `--ref_seq_input_dir`,
- read sequences through Parameter `--raw_seq_input_dir`,
- both reference and read sequences, or
- a ZIP archive through Parameter `--input_zip`.

Parameter `--graph_prefix` is used as the output prefix for the newly
created graph in this mode.

For example, a graph built with the prefix *salmonella_graph* produces files such as

- *salmonella_graph.gfa*,
- *salmonella_graph.color.bfg*, and
- *salmonella_graph.idx*.

### Build an Index for an Existing Graph

If an existing graph is provided through Parameter `--input_graph_dir` and `--graph_prefix`, the workflow creates a PLAST search index for that graph.

The directory specified through `--input_graph_dir` must contain one of either

- *&lt;prefix&gt;.gfa*, or
- *&lt;prefix&gt;.gfa.gz*, 

**and**

- *&lt;prefix&gt;.color.bfg*,

where *&lt;prefix&gt;* is the shared file name prefix specified through Parameter `--graph_prefix`.

An existing `.idx` file may be overwritten when a new index is built.

## Search

The Search mode searches an existing indexed graph with one or more query sequences.

The following parameters/inputs are required:

- Parameter `--input_graph_dir`,
- Parameter `--graph_prefix`,
- a graph file with suffix *.gfa* or *.gfa.gz*;
- a graph color file with suffix *.color.bfg*;
- a PLAST search index file with suffix *.idx*;
- Parameter `--query_fasta`.

The workflow does not build a new graph in this mode.

If one of the required graph files is missing, the workflow terminates with an error.

## Build and Search

The Build and Search mode prepares a graph and index and then performs a PLAST search.

The workflow behaves as follows:

| Available Input | Workflow Action |
|---|---|
| Reference and/or read sequence input | Build a new graph and index, then perform a search. |
| ZIP archive containing sequence input | Extract the archive, build a new graph and index, then perform a search. |
| Existing graph without sequence input | Build or rebuild the index, then perform a search. |
| Existing graph and sequence input | Build a new graph from the sequence input, then perform a search. |
| No graph and no sequence input | Terminate with an error. |

A query sequence file is required in this mode.

# Query Sequence Input

Query sequences are selected with the parameter `--query_fasta`. The input must be a [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) file containing one or more DNA query sequences.

Typical accepted file names are, e.g,

- *queries.fa*,
- *queries.fasta*, and
- *queries.fna*.

A query file may contain multiple FASTA records:

```text
>query_1
ACGTTGCAACGTTGCA

>query_2
TTGCGATCGATCGATC
```

Each FASTA record is processed as a separate PLAST query.

## Query Normalization

PLAST expects query sequences in a plain-text format with one sequence per line. FASTA headers are not passed directly to PLAST.

Before the PLAST search, the workflow normalizes the query file as follows:

1. FASTA headers are stored for later annotation of the PLAST results.
2. Multiple sequence lines belonging to one FASTA record are joined.
3. Lowercase canonical DNA bases are converted to uppercase.
4. Valid IUPAC ambiguity characters are replaced by compatible canonical DNA bases.
5. Sequences containing characters outside the accepted DNA/IUPAC alphabet are discarded.

The normalized queries are written to a separate file. Thus, the input file specified through `--query_fasta` itself is not modified.

## Accepted Canonical DNA Bases

The canonical DNA alphabet accepted by PLAST consists of:

```text
A
C
G
T
```

Lowercase versions of these bases are converted to uppercase.

For example:

```text
acgT
```

is normalized to:

```text
ACGT
```

## IUPAC Ambiguity Characters

The workflow accepts the following IUPAC ambiguity characters:

| IUPAC Character | Possible Replacement Bases |
|---|---|
| `R` | `A` or `G` |
| `Y` | `C` or `T` |
| `S` | `G` or `C` |
| `W` | `A` or `T` |
| `K` | `G` or `T` |
| `M` | `A` or `C` |
| `B` | `C`, `G`, or `T` |
| `D` | `A`, `G`, or `T` |
| `H` | `A`, `C`, or `T` |
| `V` | `A`, `C`, or `G` |
| `N` | `A`, `C`, `G`, or `T` |

For example, `R` is replaced only by `A` or `G`, because `R` represents a purine base.

The random replacement is controlled by the parameter `--random_seed`.
Using the same query file and the same value for Parameter `--random_seed` makes the normalization reproducible.

## Invalid Query Characters

Characters outside the accepted DNA/IUPAC alphabet cause the affected query sequence to be discarded.

Examples of invalid characters include:

```text
-
?
*
X
Z
1
2
```

If all query sequences are discarded, the workflow terminates with an error.

# Existing Graph Input

An existing graph is selected using Parameters `--input_graph_dir` and `--graph_prefix`.

Parameter `--input_graph_dir` specifies the bucket folder containing the graph files. Parameter `--graph_prefix` is the common prefix of all files belonging to the graph.

For example, if the selected graph directory contains:

```text
graphs/
├── myPangenome.gfa.gz
├── myPangenome.color.bfg
└── myPangenome.idx
```

then the parameter values are:

```text
--input_graph_dir = graphs/
--graph_prefix = myPangenome
```

The graph must contain exactly one graph sequence file 
*&lt;prefix&gt;.gfa* **or** *&lt;prefix&gt;.gfa.gz*, where *&lt;prefix&gt;* is the shared file name prefix of all graph files specified through Parameter `--graph_prefix`.

In addition, it must contain the graph color file *&lt;prefix&gt;.color.bfg*.

For `Search` mode, a PLAST search index file named *&lt;prefix&gt;.idx* is required.

> **Important:** PLAST requires a [Bifrost](https://github.com/pmelsted/bifrost)-compatible graph. A general GFA file
> generated by another graph tool is not necessarily suitable for PLAST.

# Sequence Directory Input

A new graph can be built from reference sequences, read sequences, or both.

The workflow provides two optional input directory parameters `--ref_seq_input_dir` and `--raw_seq_input_dir`.

At least one of these parameters must be selected when building a graph from sequence data.

## Reference Sequence Directory

The parameter `--ref_seq_input_dir` is intended for

- assembled genomes,
- contigs, and 
- reference sequences.

Files from this directory are passed to PLAST as reference input. Each input file represents one graph color.

A typical reference directory may contain:

```text
references/
├── genome_01.fasta
├── genome_02.fasta
├── genome_03.fasta.gz
└── selected_genomes.fof
```

## Read Sequence Directory

The parameter `--raw_seq_input_dir` is intended for

- raw sequencing reads, and
- other unassembled DNA sequencing data.

Files from this directory are passed to PLAST as read input and *k*-mers extracted from them will only be considered for graph construction if found at least twice.

Each input file represents one graph color.

A typical read directory may contain:

```text
reads/
├── sample_01.fastq
├── sample_02.fastq.gz
└── sample_03.fq.gz
```

## Supported Sequence File Extensions

The workflow recognizes the following sequence file extensions:

```text
.fa
.fasta
.fna
.fq
.fastq
.fa.gz
.fasta.gz
.fna.gz
.fq.gz
.fastq.gz
```

Only files directly inside the selected directory are considered. Files in subdirectories are not automatically included.

# File-of-files Input

A file-of-files file has the suffix *.fof*.

If a *.fof* file is present in a sequence input directory, the workflow uses only the sequence files listed in that file.

Each non-empty line of the *.fof* file must contain the name of one sequence file located in the same directory.

For example, the file *selected_genomes.fof* in the following directory:

```text
references/
├── selected_genomes.fof
├── genome_01.fasta
├── genome_02.fasta
├── genome_03.fasta
└── genome_04.fasta
```

may contain the content:

```text
genome_01.fasta
genome_03.fasta
genome_04.fasta
```

In this example, the workflow uses:

```text
genome_01.fasta
genome_03.fasta
genome_04.fasta
```

The file `genome_02.fasta` is ignored.

If a file listed in the `.fof` file does not exist, the workflow terminates with an error.

It is recommended to place only one `.fof` file in each sequence input
directory.

# ZIP Archive Input

Instead of selecting sequence directories separately, input data may be
provided as a ZIP archive using Parameter `--input_zip`.

The ZIP archive may contain one or both of the following directories:

```text
references/
reads/
```

At least one of these directories must be present.

A valid ZIP archive may have the following structure:

```text
plast_input.zip
├── references/
│   ├── genome_01.fasta
│   ├── genome_02.fasta
│   └── selected_genomes.fof
└── reads/
    ├── sample_01.fastq.gz
    └── sample_02.fastq.gz
```

An archive containing only reference sequences is also valid:

```text
reference_input.zip
└── references/
    ├── genome_01.fasta
    └── genome_02.fasta
```

Likewise, an archive containing only read sequences is valid, too:

```text
read_input.zip
└── reads/
    ├── sample_01.fastq.gz
    └── sample_02.fastq.gz
```

The same rules for supported sequence files and `.fof` files apply to
directories extracted from a ZIP archive.

# Search Color Set

The optional parameter `--search_color_set` restricts the PLAST search to selected graph colors.

The selected file must contain one color name per line.

A color usually corresponds to one input sequence file used during graph
construction.

For example:

```text
genome_01.fasta
genome_03.fasta
genome_07.fasta
```

When no search color set is provided, PLAST ignores colors during search.

# Parameter Overview

## Essential Parameters

| Parameter | Description |
|---|---|
| `--mode` | Selects *Build*, *Search*, or *Build and Search*. |
| `--outdir` | Bucket and folder in which workflow results are written. |
| `--query_fasta` | FASTA file containing query sequences. Required for search modes. |
| `--input_graph_dir` | Directory containing an existing graph. |
| `--graph_prefix` | Prefix identifying graph-related files. |
| `--ref_seq_input_dir` | Directory containing reference sequence files. |
| `--raw_seq_input_dir` | Directory containing read sequence files. |
| `--input_zip` | ZIP archive containing *references/* and/or *reads/*. |

## Graph Construction Parameters

| Parameter | Default | Description |
|---|---:|---|
| `--kmer_length` | 31 | Length of *k*-mers a graph is constructed or shall be constructed for. |
| `--min_length` | 23 | Length of minimizers a graph is constructed or shall be constructed for. |
| `--seed_length` | 11 | Minimal seed length used for index construction and search. |
| `--advanced_index` | *false* | Enables construction of an advanced PLAST index for faster quorum searches. |

### *k*-mer Length

The parameter `--kmer_length` defines the *k*-mer size a graph is constructed or shall be constructed for. Specification for an existing graph is optional.

A larger *k*-mer length can increase specificity but may require additional memory and may reduce sensitivity for short or highly variable sequences.

The workflow accepts values in the range from 9 to 63.

### Minimizer Length

The parameter `--min_length` defines the minimizer length a graph is constructed or shall be constructed for. Specification for an existing graph is optional.

The workflow accepts values in the range from 4 to *k* - 2, where *k* is the *k*-mer length.

### Seed Length

The parameter `--seed_length` defines the minimal exact seed length used by PLAST to initiate an alignment extension.

Allowed values are 11, 12, ..., 15.

Larger values can increase specificity but may decrease sensitivity.

### Advanced Index

The parameter `--advanced_index` enables construction of an advanced PLAST search index.

This option may improve performance for searches using a quorum restriction (cmp. Section [Quorum](#quorum)).

# Alignment and Search Parameters

| Parameter | Default | Description |
|---|---:|---|
| `--match` | 1 | Score assigned to matching bases. |
| `--mismatch` | -1 | Score assigned to mismatching bases. |
| `--gap_score` | -2 | Linear gap score used for gap opening and gap extension. |
| `--x_dropoff` | 3 | X-dropoff parameter for alignment extension. |
| `--strand` | *both* | DNA strand(s) considered during search. |
| `--max_results` | 250 | Maximum number of reported alignments per query. |
| `--quorum` | 1 | Minimum number of colors supporting an alignment. |
| `--evalue_threshold` | 10 | Maximum e-value for reported alignments. |
| `--report_colors` | *false* | Enables color coverage information in the PLAST output. |
| `--random_seed` | 42 | Seed used when replacing IUPAC ambiguity characters. |

## Gap Score

The workflow uses a linear gap score. The score may be modified using Parameter `--gap_score`.

The gap score must be negative.

## Strand Selection

The parameter `--strand` controls which DNA strand is searched.

Possible values are:

| Value | Meaning |
|---|---|
| `+` | Search only the forward strand. |
| `-` | Search only the reverse-complement strand. |
| `both` | Search both strands. |

## Maximum Number of Results

The parameter `--max_results` limits the number of reported alignments for each query sequence.

If more matching alignments are found, only the highest-scoring alignments up to the selected limit are reported.

## Quorum

The parameter `--quorum` specifies the minimum number of graph colors that must support an alignment.

A graph color generally corresponds to one sequence input file used during graph construction.

The quorum is an absolute number, not a percentage.

For example, a quorum of 10 requires a *k*-mer of the graph sequence involved in an alignment to be supported by at least ten graph colors (i.e. it appears in ten individual genomes).

## E-value Threshold

The parameter `--evalue_threshold` controls statistical filtering of reported alignments.

Lower values result in stricter filtering.

## Color Coverage Reporting

The parameter `--report_colors` enables reporting of graph color coverage for individual alignments.

This may help identify which input genomes or samples support a particular alignment.

# Valid Input Combinations

The following input combinations are supported.

| Mode | Existing Graph | Sequence Directories / ZIP Archive | Query File |
|---|---:|---:|---:|
| *Build* | Optional | Optional, but either graph or sequence input is required | Not required |
| *Search* | Required | Not required | Required |
| *Build and Search* | Optional | Optional, but either graph or sequence input is required | Required |

Examples of valid conceptual configurations:

| Intended Task | Mode | Required Inputs |
|---|---|---|
| Build a graph from references | *Build* | `--ref_seq_input_dir` |
| Build a graph from reads | *Build* | `--raw_seq_input_dir` |
| Build a graph from references and reads | *Build* | `--ref_seq_input_dir` and `--raw_seq_input_dir` |
| Build an index for an existing graph | *Build* | `--input_graph_dir` and `--graph_prefix` |
| Search an existing graph | *Search* | `--input_graph_dir`, `--graph_prefix`, and `--query_fasta` |
| Build a graph and search it | *Build and Search* | Sequence input and `--query_fasta` |
| Rebuild an index and search an existing graph | *Build and Search* | `--input_graph_dir`, `--graph_prefix`,  and `--query_fasta` |

# Invalid Input Combinations

The workflow terminates with an error in situations such as:

| Situation | Reason |
|---|---|
| Search mode without `--query_fasta` | A search requires query sequences. |
| Search mode without `--input_graph_dir` | A search requires an existing graph. |
| Existing graph without `--graph_prefix` | The graph files cannot be identified. |
| Existing graph has no *.gfa* or *.gfa.gz* file | Graph sequence information is missing. |
| Existing graph has no *.color.bfg* file | Graph color information is missing. |
| Search mode without *.idx* file | PLAST searching requires an index. |
| Build mode without graph input, sequence directories, or ZIP archive | No graph can be built or indexed. |
| Build and Search mode without graph input, sequence directories, or ZIP archive | No graph can be created or selected. |
| Build and Search mode without `--query_fasta` | The final search cannot be performed. |
| No valid query remains after normalization | PLAST cannot perform a search without valid query sequences. |
| Invalid `--kmer_length`, `--min_length`, `--seed_length`, `--gap_score`, or `--quorum` values | Parameter constraints are violated. |

# Output

The parameter `--outdir` defines where workflow results are written.

In the CloWM parameter form, select

1. an output bucket, and
2. a folder inside that bucket.

The workflow writes all final result files to this selected location.

After the workflow has finished:

1. Navigate to **Files > My Data Buckets**.
2. Select the output bucket.
3. Open the output folder as specified through `--outdir`.
4. Download the desired result files.

Refer to Tab `Output` for more information about generated outputs.

# Monitoring Workflow Execution

After launching the workflow, CloWM opens the workflow execution report.

The report provides information about

- the current execution status,
- start time and end time,
- running and completed workflow processes,
- CPU and memory usage,
- input/output activity, and
- process-specific logs after workflow completion.

If a workflow fails, inspect the displayed logs and the Nextflow log files
available in the execution report.

# Local Execution

The workflow can also be executed locally with Nextflow and Docker. Please refer to [PLAST's code repository on Github](https://github.com/gi-bielefeld/plast) for the source code!
