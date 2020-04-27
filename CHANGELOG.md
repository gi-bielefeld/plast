# Changelog

* v0.1.0 (04-27-2020)
  - Pure quorum checks have been made faster now by two changes: Before a matrix
    is iterated over to check fulfillment of quorum *m*, the maximum color id 
    *c* being present in that matrix is queried. If already $`c <= m`$ holds, 
    the check is aborted. 
    Additionally, missing colors are now counted next to present ones as well 
    which allows to terminate checks also if too many colors have been missed 
    already.