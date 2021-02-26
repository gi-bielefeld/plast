#ifndef IO_HPP
#define IO_HPP

#include "Search.h"

#define MIN_PARAM_NB 4
#define OPTIONS "i:s:w:k:g:S:R:t:q:c:Q:m:M:X:e:n:o:d:l:L:C:T:aru"
#define BUILD_COMMAND "Build"
#define SEARCH_COMMAND "Search"
#define RUNTIME_FLAG_DEFAULT false
#define PREPROS_FLAG_DEFAULT -1
#define GFA_FILE_ENDING ".gfa"
#define COLOR_FILE_ENDING ".bfg_colors"
#define INDEX_FILE_ENDING ".idx"
#define SEARCH_STRAND_DEFAULT Both
#define PLUS_STRAND '+'
#define MINUS_STRAND '-'
#define BASES_PER_LINE 60
#define GAP '-'
#define REPORT_COLORS_FLAG_DEFAULT false

//This function reads in a file in which colors are stored the search will be based on
const list<pair<string, size_t>> loadSearchColors(const char* filename, uint32_t& nbCols);

//This function parses the program parameters. Returns false if given arguments are not valid
const bool parseArgs(int& nb_args, char** argList, int16_t& prepros, string& filePref, int32_t& s, int32_t& k, int32_t& g, CCDBG_Build_opt &gOpt, int32_t& t, string& qFile, string& c, uint32_t& m, SrchStrd& strd, bool& r, uint16_t &mscore, int16_t &mmscore, int16_t& X, int32_t &gOpen, int32_t &gExt, uint16_t &nRes, double &lambda, double &lambdaG, double &C, double &Cgap, double &eValLim, bool &isSim, bool &advIdx);

//This function prints usage infos
inline void dispHelp(){
	cerr << "PLAST [COMMAND] [COMMAND_PARAMETERS]" << endl << endl;
	cerr << "[COMMAND]:" << endl << endl;
	cerr << "   Build    Build index (and graph)" << endl;
	cerr << "   Search   Search inside an indexed graph" << endl << endl;
	cerr << "[COMMAND_PARAMETERS]: Build" << endl << endl;
	cerr << "   >Mandatory with required argument:" << endl << endl;
	cerr << "   -i   --graph-prefix   File prefix of a graph which already exists or has to be built" << endl << endl;
	cerr << "   >Optional with required argument:" << endl << endl;
	cerr << "   -w   --seed-length   Minimal seed length an index is built for (default is 11)" << endl;
	cerr << "   -k   --kmer-length   Length of k-mers in a newly built graph (default is 31)" << endl;
	cerr << "   -g   --min-length    Length of minimizers in a newly built graph (default is 23)" << endl;
	cerr << "   -S   --input-seqs    Names of raw input sequence file(s) (FASTA/FASTQ) to build a new graph from (all sequences from the same file will share a color in the graph)" << endl;
	cerr << "   -R   --input-refs    Names of reference input sequence file(s) (FASTA/FASTQ) to build a new graph from (all sequences from the same file will share a color in the graph)" << endl;
	cerr << "   -t   --threads       Number of threads to be used for graph construction and loading (default is 1)" << endl << endl;
	cerr << "   >Optional without argument:" << endl << endl;
	cerr << "   -a   --advanced-index   Construct advanced index for faster quorum searches" << endl << endl;
	cerr << "[COMMAND_PARAMETERS]: Search" << endl << endl;
	cerr << "   >Mandatory with required argument:" << endl << endl;
	cerr << "   -q   --query          Query sequence file (1 query per line)" << endl;
	cerr << "   -i   --graph-prefix   Graph, color and index file prefix" << endl << endl;
	cerr << "   >Optional with required argument:" << endl << endl;
	cerr << "   -w   --seed-length    Minimal seed length (default is 11)" << endl;
	cerr << "   -s   --search-set     Search color set file with one color name per line to consider during the search (default is all colors)" << endl;
	cerr << "   -Q   --quorum         Absolute quorum (default is 1 color)" << endl;
	cerr << "   -m   --mismatch       Mismatch score used (default is -1)" << endl;
	cerr << "   -M   --match          Match score used (default is 1)" << endl;
	cerr << "   -X   --X-dropoff      X-dropoff value" << endl;
	cerr << "   -d   --gap-open       Gap open score used for gapped alignment (default is -2)" << endl;
	cerr << "   -e   --gap-extension  Gap extension score used for gapped alignment (default is -2)" << endl;
	cerr << "   -n   --max-results    Maximum number of alignments to be outputted (default is " << DEFAULT_NB_RES << ")" << endl;
	cerr << "   -o   --strand         DNA strands to consider during search. Can be '+','-' (default is both)" << endl;
	cerr << "   -l   --lambda         Statistical value lambda for ungapped extension" << endl;
	cerr << "   -L   --lambda-gap     Statistical value lambda for gapped extension" << endl;
	cerr << "   -c   --stat-C         Statistical value C for ungapped extension" << endl;
	cerr << "   -C   --stat-C-gap     Statistical value C for gapped extension" << endl;
	cerr << "   -T   --e-val-thres    Expectation value threshold for a hit to be considered (default is 10)" << endl << endl;
	cerr << "   >Optional without argument:" << endl << endl;
	cerr << "   -r   --report-colors   Enable alignment color coverage output" << endl;
}

//This function loads query sequences from a file and stores it in a list
void loadQueries(const string &filename, vector<string> &qList);

//This function outputs a hit's alignment
void repAlgn(const Hit *res);

//This function outputs the color sets of a given result
void outpColSets(ColoredCDBG<UnitigInfo> &cdbg, const Hit *res);

//This function iterates over an UnitigColors object and returns the set of all color names
inline const vector<string> formColSet(ColoredCDBG<UnitigInfo> &cdbg, const UnitigColorMap<UnitigInfo> &uni){
	string color;
	//Initialize color set
	vector<string> set = vector<string>();
	//Get UnitigColors object
	const UnitigColors uniCols = *uni.getData()->getUnitigColors(uni);

	//Iterate over colors
	for(UnitigColors::const_iterator c = uniCols.begin(uni); c != uniCols.end(); ++c){
		//Get current ID's color name
		color = cdbg.getColorName(c.getColorID());

		//If the color set is empty we have nothing to compare against
		if(set.empty()){
			set.push_back(color);
		} else if(set.back() != color){//Check if we have just seen the current color already
			set.push_back(color);
		}
	}

	return set;
}

#endif