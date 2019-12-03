# PLAST

### Pangenome Local Alignment Search Tool

* Local alignment search between DNA query and graphical pangenome 
* Assembled genomes or reads as input
* Significance estimation based on alignment statistic
* Customized focus on parts of the graph
* Color coverage output per alignment on demand

## Table of Contents

* [Requirements](https://gitlab.ub.uni-bielefeld.de/gi/plast#requirements)
* [Compilation](https://gitlab.ub.uni-bielefeld.de/gi/plast#compilation)
* [Usage](https://gitlab.ub.uni-bielefeld.de/gi/plast#usage)
* [Test data](https://gitlab.ub.uni-bielefeld.de/gi/plast#testdata)
* [FAQ](https://gitlab.ub.uni-bielefeld.de/gi/plast#faq)
* [Contact](https://gitlab.ub.uni-bielefeld.de/gi/plast#contact)
* [License](https://gitlab.ub.uni-bielefeld.de/gi/plast#license)

## Requirements

PLAST builds a **compacted, colored de Bruijn graph** from given input genomes using the API of [Bifrost](https://github.com/pmelsted/bifrost). Apart from the requirements of Bifrost (c++ and cmake), there are no further strict dependencies.

Parameter estimations for alignment statistics can be done using a workflow that
requires a running version of 
[snakemake](https://snakemake.readthedocs.io/en/stable/) and the additional
packages [matplotlib](https://matplotlib.org) and
[scipy](https://www.scipy.org).

## Compilation

```
cd <plast_directory>/src
make
```

By default, the installation creates:
* a binary (*PLAST*)

You may want to make the binary (*PLAST*) accessible via your (*PATH*) variable.

Please note the installation instructions regarding the default maximum *k*-mer size of Bifrost from its README.
If your Bifrost libraries have been compiled for 64 bit, change the PLAST 
makefile accordingly (add `-DMAX_KMER_SIZE=64` to CFLAGS).

If during the compilation, the bifrost library files are not found, make sure that the corresponding folder is found as include path by the C++ compiler. You may have to add
`-I/usr/local/include` (with the corresponding folder) to the CFLAGS in the makefile.

## Usage:

```
PLAST
```

displays the command line interface:
```
PLAST [COMMAND] [COMMAND_PARAMETERS]

[COMMAND]:

   Build    Build index (and graph)
   Search   Search inside an indexed graph

[COMMAND_PARAMETERS]: Build

   >Mandatory with required argument:

   -i   --graph-prefix   File prefix of a graph which already exists or has to be built

   >Optional with required argument:

   -w   --seed-length   Minimal seed length an index is built for (default is 11)
   -k   --kmer-length   Length of k-mers in a newly built graph (default is 31)
   -g   --min-length    Length of minimizers in a newly built graph (default is 23)
   -S   --input-seqs    Names of raw input sequence file(s) (FASTA/FASTQ) to build a new graph from (all sequences of a file will share a color in the graph)
   -R   --input-refs    Names of reference input sequence file(s) (FASTA/FASTQ) to build a new graph from (all sequences of a file will share a color in the graph)
   -t   --threads       Number of threads to be used for graph construction (default is 1)

[COMMAND_PARAMETERS]: Search

   >Mandatory with required argument:

   -q   --query          Query sequence file (1 query per line)
   -i   --graph-prefix   Graph, color and index file prefix

   >Optional with required argument:

   -w   --seed-length    Minimal seed length (default is 11)
   -s   --search-set     Search color set file with one color name per line to consider during the search (default is all colors)
   -m   --quorum         Quorum (absolute value)
   -X   --X-dropoff      X-dropoff value
   -n   --max-results    Maximum number of alignments to be outputted (default is 250)
   -d   --strand         DNA strands to consider during search. Can be '+','-' (default is both)
   -l   --lambda         Statistical value lambda for ungapped extension
   -L   --lambda-gap     Statistical value lambda for gapped extension
   -c   --stat-C         Statistical value C for ungapped extension
   -C   --stat-C-gap     Statistical value C for gapped extension
   -e   --e-value        Expectation value threshold for a hit to be considered (default is 10)

   >Optional without argument:

   -r   --report-colors   Enable alignment color coverage output
```

### Examples

1. **Graph building and indexing**

   If a graph has to be built and indexed from scratch, a call might look like

   ```
   PLAST Build -i someGraph -R genomeAssemblyA.fasta genomeAssemblyB.fasta -S genomeReadsC.fasta genomeReadsD.fasta -t 4
   ```

   The colored de-Bruijn graph is built with Bifrost using 4 threads (`-t 4`) from assemblies *genomeAssemblyA.fasta* and *genomeAssemblyB.fasta* (`-R`) and read datasets *genomeReadsC.fasta* and *genomeReadsD.fasta* (`-S`). Graph, color information and index are saved as *someGraph.{gfa,bfg_colors,idx}* (`-i someGraph`).

2. **Indexing of existing graph**

   If a colored de-Bruijn graph has already been built using Bifrost, it can just be indexed by

   ```
   PLAST Build -i anotherGraph -k 25 -g 18 -w 12
   ```

   The index is built using the graph saved as *anotherGraph.{gfa,bfg_colors}* 
   (`-i anotherGraph`) which has been built with a k-mer length of 25 (`-k 25`) 
   and minimizer length of 18 (`-g 18`) for a minimal seed length of 12 
   (`-w 12`). The index is saved as *anotherGraph.idx*.
   
   **ATTENTION:** Be aware that for using values of k > 31 your Bifrost libraries 
   have to be compiled for 64 bit and `-DMAX_KMER_SIZE=64` has to be added to 
   the CFLAGS in PLAST's makefile before compilation!

3. **Searching**

   A basic search may look like this:

   ```
   PLAST Search -i someGraph -q myQueries.q
   ```

   We want to search three queries stored in *myQueries.q* (one query per line) within the graph *someGraph.{gfa,bfg_colors,idx}* (see Example 1).

   Assuming we are only interested in alignments to genomic sequences of either 
   *genomeAssemblyA.fasta*, *genomeAssemblyB.fasta* or *genomeReadsC.fasta*. In 
   the context of a colored de-Bruijn graph, sequences from these input files 
   each have their own color in the graph. We have to specify the colors in a 
   search set file (`-s mySearchSet.txt`):

   ```
   PLAST Search -i someGraph -q myQueries.q -s mySearchSet.txt
   ```

   File *mySearchSet* consists of one line per color. A color's name is identical to the corresponding input file's name. In our case, it may look like:

   ```
   genomeAssemblyA.fasta
   genomeAssemblyB.fasta
   genomeReadsC.fasta
   ```

   If we additionally want to require that sequences of at least two colors 
   should support each alignment that we want to find, we have to specify this 
   by a quorum (`-m 2`):

   ```
   PLAST Search -i someGraph -q myQueries.q -s mySearchSet.txt -m 2
   ```

4. **Alignment statistic parameter estimation**
   
   PLAST uses alignment statistics in order to filter out alignments not representing sequence homology between query and graph sequence and to calculate an
   e-value for each result. Default parameters for our alignment statistic are chosen carefully and allow basic searches in any pangenome graph. However,
   results with an e-value close to the significance threshold have to be handled with care.

   Sound statistical parameters can be estimated by running a 
   [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow provided in
   the *simulation* directory. Alternatively, a Bash script is provided and may
   be used for the same purpose as well.
   
   For running the snakemake workflow, proceed as follows:
   
   1. Change to the directory *simulation*:
   
      ```
      cd simulation
      ```

   2. Before executing the workflow, it requires the location of a graph
      parameters have to be
      estimated for. Add the graph's path to the configuration file
      (*config.yaml*). You may need to change graph properties (*k*-mer and/or
      minimizer length) or the minimal seed length your graph was indexed for if
      they are different from default.


      *config.yaml*:
      ```
      ...
      #Path to the pangenome graph (please insert here!)
      gPathPref: "/path/to/myGraph.gfa"
      
      #Graph and index properties (please change if non-default!)

      #k-mer length used for graph
      k: 31
      #Minimizer length used for graph
      miniLen: 23
      #Minimal seed length used for index
      seedLen: 11
      ...
      ```
   
   3. Run the workflow by typing `snakemake`. Depending on your graph, simulations might take a while. Use the option `--cores` to run simulations on many
      cores in parallel.
   
      ```
      snakemake --cores <Nb_cores>
      ```

   4. Simulation results for gapped and ungapped alignment parameters can be found in *results/parameters.txt*.
    
   Running the workflow is the recommened way for parameter estimation. However,
   using a simple Bash script is possible as well.

   For running the Bash script, proceed as follows:

   1. Change to the directory *simulation*:
   
      ```
      cd simulation
      ```

   2. Execute the Bash script by providing the graph path as first argument.
      Non-default values for *k*-mer length (`-k`), minimizer length (`-g`) and
      minimal seed length (`-w`) may follow using the appropriate flag.

      ```
      ./estParams.sh /path/to/myOtherGraph.gfa -k <my_k> -g <my_minimizer_length> -w <my_min_seed_length>
      ```
      
   3. Simulation results for gapped and ungapped alignment parameters can be
      found in *results/parameters.txt*.

## Test data

Test data is provided in the directory *testdata*.

1. **Pangenome simulation**

   The directory contains a parameter file to be used for pangenome simulation 
   with [ALF](http://alfsim.org/#index). Additionally, a script is provided to 
   transfer sequences outputted by ALF into contiguous genome sequences.
   
3. **PARA C**

   All assembly barcodes of the Para C dataset published [here](https://genome.cshlp.org/content/28/9/1395.short) are listed in *ParaCcomplete220.txt*. In combination with an access token, they can be downloaded from 
   [EnteroBase](https://enterobase.warwick.ac.uk) using the the provided download script:
   
   ```
   ./download_assemblies.sh ParaCcomplete220.txt <token>
   ```
   
   How to get a token is described [here](https://bitbucket.org/enterobase/enterobase-web/wiki/Getting%20started%20with%20Enterobase%20API).
   
   File names of assemblies subsets of sizes 12 and 75 can be found in *ParaCsubset12.txt* and *ParaCsubset75.txt*.
   
4. **Salmonella pangenome from EnteroBase**

   A pangenome of 5000 *Salmonella typhimurium* assemblies from [EnteroBase](https://enterobase.warwick.ac.uk). A list of all 5000 assembly
   barcodes can be used to download the data from EnteroBase. For this, an access token is required. How to get a token is described
   [here](https://bitbucket.org/enterobase/enterobase-web/wiki/Getting%20started%20with%20Enterobase%20API).
   
   Use the provided script for downloading the assemblies:
   
   ```
   ./download_assemblies.sh chosen5000.txt <token>
   ```

## FAQ

We recommend to have a look at the [FAQs of Bifrost](https://github.com/pmelsted/bifrost#faq).


## Contact

For any question, feedback or problem, please feel free to file an issue on this Git repository or contact the developers and we will get back to you as soon as possible.

## License

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)

* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)

* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)

* The kseq library is copyrighted by Heng Li and released
  under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)

* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)

* Bifrost is BSD-2 licensed (https://github.com/pmelsted/bifrost)

* PLAST is GNU GPLv3 licensed [LICENSE](https://gitlab.ub.uni-bielefeld.de/gi/plast/blob/master/LICENSE)