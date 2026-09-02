# PLAST <img src="plast-logo_100.png" style="border:0;" alt="PLAST icon" align="right"/>
### Pangenome Local Alignment Search Tool

PLAST is a tool for **fast local alignment search between a nucleotide sequence query and a pangenome graph** consisting of thousands of individual genomes. 
As input, PLAST accepts prebuilt pangenome graphs as well as sets of individual genome sequences which are then used to construct a pangenome graph in a preprocessing step.
Genome sequences to be used for graph construction may be **completely assembled genome sequences, sets of contigs, raw read data**, or a mixture of such. 
Local alignments between query and all sequences represented by the pangenome graph are searched based on a seed-and-extend approach. 
**Searches may be customized** by restricting them to consider only sequences of specific genomes or sequences which occur in at least a certain number of individual genomes.
Alignments are **outputted along with a significance value** and optionally also with **information about which individual genomes** from the graph are involved in the alignment.

Further details about the method of PLAST may be found in its [publication](https://doi.org/10.1093/bioinformatics/btab077).

## Running PLAST Locally

PLAST may also be downloaded and installed locally from our [git repository](https://github.com/gi-bielefeld/plast).

## Contact

For any question, feedback or problem, please feel free to file an issue on [Github](https://github.com/gi-bielefeld/plast) or [contact the developers](mailto:pangenomics-service@cebitec.uni-bielefeld.de) and we will get back to you as soon as possible.

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
