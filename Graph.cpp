#include "Graph.h"

//This function builds a graph
void genGraph(ColoredCDBG<seedlist> &cdbg, CCDBG_Build_opt &cdbgOpt){
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

//This function assigns color ids from the graph to color names in the color set
void mapColorIds(list<pair<string, size_t>>& colorSet, const ColoredCDBG<seedlist>& g){
	//Get number of colors in the graph
	size_t totCol = g.getNbColors();

	//Iterate over color ids of the graph
	for(size_t i = 0; i < totCol; ++i){
		//Iterate over color set colors
		for(list<pair<string, size_t>>::iterator j = colorSet.begin(); j != colorSet.end(); ++j){
			//Check whether the current color id belongs to the current color name and save the info
			if(j->first == g.getColorName(i)) j->second = i;
		}
	}
}