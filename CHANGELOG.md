# Changelog

* v0.1.0 (04-27-2020)
  - Pure quorum checks have been made faster now by two changes: Before a matrix
    is iterated over to check fulfillment of quorum $`m`$, the maximum color id 
    $`c`$ being present in that matrix is queried. If already $`c <= m`$, 
    the check is aborted. 
    Additionally, missing colors are now counted next to present ones as well 
    which allows to terminate checks also if too many colors have been missed 
    already.

* v0.2.0 (03-31-2021)
  - Hits with invalid offsets are now moved to successive unitigs and kept even 
    though they do not have a right extension path.
  - The unitig on which a left, gapped extension ends is now saved for each hit 
    to improve the result filtration after gapped extension.
  - An input parameter check has been added to ensure that a given minial seed 
    length is always positive.
  - Seeds that could not be extended without gaps are no longer discarded if 
    they are longer than the minimal seed length.
  - Quorum precalculation is now possible during index building. Information is
    saved as part of the index
  - Introduced a new index format to enable storage of quorum information.
  - Enabled multi-threaded graph loading
  - Closed some smaller memory leaks.
  - Disabled verbose mode for graph reading and removed some logging info.
  - Some refactoring.
