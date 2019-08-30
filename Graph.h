#ifndef GRAPH_HPP
#define GRAPH_HPP

//Default k values for the compacted CDBG
#define DEFAULT_GRAPH_K 31//9
//Default minimizer length for the compacted CDBG
#define DEFAULT_GRAPH_G 23//4
//Default number of threads used for graph construction
#define DEFAULT_NB_THREADS 1

//This function builds a graph
void genGraph(ColoredCDBG<seedlist> &cdbg, vector<string> &seqs);

//This function assigns color ids from the graph to color names in the color set
void mapColorIds(list<pair<string, size_t>>& colorSet, const ColoredCDBG<seedlist>& g);

#endif