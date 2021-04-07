#include "Graph.h"

//This function builds a graph
void genGraph(ColoredCDBG<UnitigInfo> &cdbg, CCDBG_Build_opt &cdbgOpt){
	//Setting build options

	//Print information messages during execution if true. Default is false.
	cdbgOpt.verbose = true;
	//Prefix for the name of the file to which the graph must be written. Mandatory parameter.
	cdbgOpt.prefixFilenameOut = "";
	//Build the graph
	cdbg.build(cdbgOpt);
	cdbg.simplify(cdbgOpt.deleteIsolated, cdbgOpt.clipTips, cdbgOpt.verbose);
	//Map colors on the graph
	cdbg.buildColors(cdbgOpt);
}