#include "Graph.h"

//This function builds a graph
void genGraph(ColoredCDBG<seedlist> &cdbg, CCDBG_Build_opt &cdbgOpt){
	//Setting build options
	//cdbgOpt.outputColors = true;//Boolean indicating if color sets must be written to disk (true) or not (false). Default is true.

	//cdbgOpt.outputGFA = false;//Boolean indicating if the graph is written to a GFA file (true) or if the unitigs are written to a FASTA file (false). Default is true. NOTE: Deprecated and, thus, always true
	cdbgOpt.verbose = true;//Print information messages during execution if true. Default is false.
	//cdbgOpt.useMercyKmers = false;//Keep in the graph low coverage k-mers connecting tips of the graph. Default is false.

	//Testing
	// cout << "cdbgOpt.k:" << cdbgOpt.k << endl;
	// cout << "cdbgOpt.g:" << cdbgOpt.g << endl;

	// cdbgOpt.k = cdbg.getK();//Length of k-mers (not used by ColoredCDBG::build). Default is 31.
	// cdbgOpt.g = 46;//Length of g-mers, the minimizers, such that g < k (not used by ColoredCDBG::build).

	//Testing
	// cout << "cdbgOpt.k:" << cdbgOpt.k << endl;
	// cout << "cdbgOpt.g:" << cdbgOpt.g << endl;

	cdbgOpt.prefixFilenameOut = "";//Prefix for the name of the file to which the graph must be written. Mandatory parameter.
	
	// cdbgOpt.filename_seq_in = seqs;//vector of strings, each string is the name of a FASTA/FASTQ/GFA to file use for the graph construction. Mandatory parameter.
	//cdbgOpt.filename_seq_in.push_back("ERR431464_1.fasta");
	//cdbgOpt.filename_seq_in.push_back("ERR431464_2.fasta");
	//cdbgOpt.filename_seq_in.push_back("Test11_color2.fa");

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