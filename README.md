# PLAST

### Pangenome Local Alignment Search Tool

* Local alignment search between DNA query and graphical pangenome 
* Assembled genomes or reads as input
* Significance estimation based on alignment statistic
* Customized focus on parts of the graph
* Color coverage output per alignment on demand

### Publication

Schulz, T., Wittler, R., Rahmann, S., Hach, F., Stoye, J.: [Detecting High Scoring Local Alignments in Pangenome Graphs.](https://doi.org/10.1093/bioinformatics/btab077). Bioinformatics. (2021)

## Table of Contents

* [Requirements](https://github.com/gi-bielefeld/plast#requirements)
* [Compilation](https://github.com/gi-bielefeld/plast#compilation)
* [Usage](https://github.com/gi-bielefeld/plast#usage)
* [Test data](https://github.com/gi-bielefeld/plast#test-data)
* [Tool comparison](https://github.com/gi-bielefeld/plast#tool-comparison)
* [FAQ](https://github.com/gi-bielefeld/plast#faq)
* [Contact](https://github.com/gi-bielefeld/plast#contact)
* [Licenses](https://github.com/gi-bielefeld/plast#licenses)
* [Privacy](https://github.com/gi-bielefeld/plast#privacy)

## Requirements

PLAST builds a **compacted, colored de Bruijn graph** from given input genomes using the API of [Bifrost](https://github.com/pmelsted/bifrost). Apart from the requirements of Bifrost (c++ and cmake), there are no further strict dependencies.

Parameter estimations for alignment statistics can be done using a workflow that
requires [snakemake](https://snakemake.readthedocs.io/en/stable/) (version 
5.20.1 or higher) and the Python packages [matplotlib](https://matplotlib.org) 
and [scipy](https://www.scipy.org).

A provided tool comparison workflow additionally requires the package 
[Biopython](https://biopython.org) to be installed on your system.

## Compilation

```
cd <plast_directory>/src
make
```

By default, the installation creates:
* a binary (*PLAST*)

You may want to make the binary accessible via your *PATH* variable.

Please note the installation instructions regarding the default maximum *k*-mer size of Bifrost from its README.
E.g., if your Bifrost libraries have been compiled to support a *k*-mer length of up to 63, change the PLAST 
makefile accordingly (add `-DMAX_KMER_SIZE=64` to CFLAGS).

If during the compilation, the bifrost library files are not found, make sure that the corresponding folder is found as include path by the C++ compiler. You may have to add
`-I/usr/local/include` (with the corresponding folder) to CFLAGS in the makefile.

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
   -S   --input-seqs    Names of raw input sequence file(s) (FASTA/FASTQ) to build a new graph from (all sequences from the same file will share a color in the graph)
   -R   --input-refs    Names of reference input sequence file(s) (FASTA/FASTQ) to build a new graph from (all sequences from the same file will share a color in the graph)
   -t   --threads       Number of threads to be used for graph construction and loading (default is 1)

   >Optional without argument:

   -a   --advanced-index   Construct advanced index for faster quorum searches

[COMMAND_PARAMETERS]: Search

   >Mandatory with required argument:

   -q   --query          Query sequence file (1 query per line)
   -i   --graph-prefix   Graph, color and index file prefix

   >Optional with required argument:

   -w   --seed-length    Minimal seed length (default is 11)
   -s   --search-set     Search color set file with one color name per line to consider during the search (default is all colors)
   -Q   --quorum         Absolute quorum (default is 1 color)
   -m   --mismatch       Mismatch score used (default is -1)
   -M   --match          Match score used (default is 1)
   -X   --X-dropoff      X-dropoff value
   -d   --gap-open       Gap open score used for gapped alignment (default is -2)
   -e   --gap-extension  Gap extension score used for gapped alignment (default is -2)
   -n   --max-results    Maximum number of alignments to be outputted (default is 250)
   -o   --strand         DNA strands to consider during search. Can be '+','-' (default is both)
   -l   --lambda         Statistical value lambda for ungapped extension
   -L   --lambda-gap     Statistical value lambda for gapped extension
   -c   --stat-C         Statistical value C for ungapped extension
   -C   --stat-C-gap     Statistical value C for gapped extension
   -T   --e-val-thres    Expectation value threshold for a hit to be considered (default is 10)

   >Optional without argument:

   -r   --report-colors   Enable alignment color coverage output
```

### Examples

1. **Get some testing data**

   You may want to get some toy data first to try out the program.

   ```
   wget 'https://hgdownload.cse.ucsc.edu/goldenPath/eboVir3/bigZips/160sequences.tar.gz'
   tar xvzf 160sequences.tar.gz
   ```

   These commands let you download and unpack a small dataset of 160 individual samples of different Ebola virus subspecies.

1. **Graph building and indexing**

   In order to build a graph from this input data and indexing it, a call might look like

   ```
   PLAST Build -i ebolaPangenome -R *.fa -t 4
   ```

   The above command builds a colored de-Bruijn graph via the Bifrost library using 4 threads (`-t 4`) from all Ebola sequences. Graph, color information and index are saved as *ebolaPangenome.{gfa,bfg_colors,idx}* (`-i ebolaPangenome`). Note that our whole dataset consists of already assembled data. Thus, they are passed to the program using `-R`. Assuming we would deal with raw read data here, you would rather want to use `-S` to allow some quality filtering during graph construction. Even a combination of these two parameters is possible if dealing with both kinds of data.

   The graph is build using default _k_-mer size of 31. This may be changed using option `-k <KMER-LENGTH>`.

   **ATTENTION:** Be aware that for using values of _k_ > 31 your Bifrost library 
   has to be compiled using the option `-DMAX_KMER_SIZE=x`, where *x* has to be a multiple of 32. The option has to be added to 
   CFLAGS in PLAST's makefile before compilation, too!

2. **Indexing of an existing graph**

   The above command has built an index using the default minimal seed length of 11. Assuming we would want to use a different minimal seed length from now on, a new index may also be calculated for an already existing Bifrost graph. The command

   ```
   PLAST Build -i ebolaPangenome -w 12
   ```

   recalculates an index for our Ebola graph (`-i ebolaPangenome`) using a minimal seed length of 12 (`-w 12`). It automatically overwrites the old index file.

3. **Searching**

   Inside the directory `testdata`, you may find a file containing two _unknown_ sequences (one sequence per line). 
   A basic search for these sequences inside the Ebola graph would look like

   ```
   PLAST Search -i ebolaPangenome -q unknownQueries.q -w 12
   ```

   Alignment results are outputted on the command line seperately for each query. They reveal a perfect match for each query and various suboptimal ones. 
   
   Repeating the search a second time using parameter `-r` allows us to get information about the individual samples from the graph involved in each alignment. Both perfect matches (first result for each query) are completely covered by a sample of Marburg virus (*NC_024781v1*) and partially by a second sample (*NC_001608v3*) as well.

   Assuming we are particularly interested in how alignments look for NC_001608v3, we could specify a *search set* for our search. It allows to focus alignment searches inside the graph exclusively to sequences from certain samples. A search set may be passed to the algorithm using a text file with one *sample id* per line. A sample id here means the sample's file path used for graph building.
   A search set for our current use case may be found below the *testdata* directory, too.

   The command

   ```
   PLAST Search -i ebolaPangenome -q unknownQueries.q -w 12 -r -s searchSet.txt
   ```

   now reveals alignments completely covered by NC_001608v3 but lower scores.

   Otherwise, we could also be interested only in alignments that involve at least a certain amount of all samples. This may be realized using a *quorum*.
   
   The command 

   ```
   PLAST Search -i ebolaPangenome -q unknownQueries.q -w 12 -Q 80
   ```

   allows to find all alignments supported by at least 80 (50%) of all samples of our Ebola graph.

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
      (*config.yaml*). You may need to change graph properties (_k_-mer and/or
      minimizer length) or the minimal seed length your graph was indexed for if
      they are different from default.


      *config.yaml*:
      ```
      ...
      #Path to the pangenome graph (please insert here!)
      gPathPref: "/path/to/graph.gfa"
      
      #Graph and index properties (please change if non-default!)

      #k-mer length used for graph
      k: 31
      #Minimizer length used for graph
      miniLen: 23
      #Minimal seed length used for index
      seedLen: 11
      ...
      ```
   
   3. Run the workflow by typing `snakemake`. Depending on your graph, simulations might take a while. Use the option `--cores` to specify the number of cores the simulation shall be performed on in parallel.
   
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
      ./estParams.sh /path/to/graph.gfa -k 31 -g 23 -w 11
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

   All assembly barcodes of the Para C dataset published [here](https://genome.cshlp.org/content/28/9/1395.short) are listed in *ParaCcomplete220.txt*. Data are provided from 
   [EnteroBase](https://enterobase.warwick.ac.uk). The exact sequences used in our publication can also be downloaded [here](https://gitlab.ub.uni-bielefeld.de/gi/parac).
   
   File names of assembly subsets of sizes 12 and 75 can be found in *ParaCsubset12.txt* and *ParaCsubset75.txt*.
   
4. **Salmonella pangenome from EnteroBase**

   A pangenome of 5000 *Salmonella typhimurium* assemblies from [EnteroBase](https://enterobase.warwick.ac.uk) can be downloaded [here](https://gitlab.ub.uni-bielefeld.de/gi/typhimurium). A list of the exact 5000 assemblies we used in our publication can be found in *chosen5000.txt*.

## Tool comparison

A run time and memory usage comparison is documented as a 
[snakemake](https://snakemake.readthedocs.io/en/stable/) workflow in the 
directory 
*comparison*.

For execution of the workflow, proceed as follows:

* Install all programs to be tested on your system.

* Get some testing data.

* Provide the workflow with the locations of your testing data and program 
  binaries by editing the configuration file *comparison/config.yaml*:

  ```
  ...
  # PLEASE ADJUST THE FOLLOWING PARAMETERS --------------------------------------

  #The place where the data is stored
  dataDir: "path/to/my_testing_data_dir"

  #The program binaries that shall be used
  blastnbin: "path/to/blastn_binary"
  makeblastdbBin: "path/to/makeblastdb_binary"
  mmseqs2Bin: "path/to/mmseqs2_binary"
  ublastBin: "path/to/usearch_binary"
  #Only an absolute path will work here!
  blatBin: "/absolute/path/to/blat_binary"
  blat_faToTwoBit: "/absolute/path/to/faToTwoBit_binary"
  ...
  ```
  
* Change into directory *comparison* and run the workflow.

  ```
  cd comparison
  snakemake
  ```

## FAQ

We recommend to have a look at the [FAQs of Bifrost](https://github.com/pmelsted/bifrost#faq).


## Contact

For any question, feedback or problem, please feel free to file an issue on [Github](https://github.com/gi-bielefeld/plast) or [contact](mailto:pangenomics-service@cebitec.uni-bielefeld.de) the developers and we will get back to you as soon as possible.

PLAST is provided as a service of the [German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/). We would appriciate if you would participate in the evaluation of PLAST by completing this [very short survey](https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=plast).

## Licenses

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)

* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)

* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)

* The kseq library is copyrighted by Heng Li and released
  under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)

* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)

* Bifrost is BSD-2 licensed (https://github.com/pmelsted/bifrost)

* PLAST is GNU GPLv3 licensed [LICENSE](https://github.com/gi-bielefeld/plast/blob/master/LICENSE)

## Privacy

We use the open source software Matomo for web analysis in order to collect anonymized usage statistics for this repository. Please refer to our [Privacy Notice](/PrivacyNotice.pdf) for details.

<!-- Matomo Image Tracker-->
<img referrerpolicy="no-referrer-when-downgrade" src="https://piwik.cebitec.uni-bielefeld.de/matomo.php?idsite=28&amp;rec=1" style="border:0" alt="" />
<!-- End Matomo -->
