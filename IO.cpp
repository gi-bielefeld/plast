#include "IO.h"
#include <getopt.h>

//This function parses the program parameters. Returns false if given arguments are not valid.
const bool parseArgs(int& nb_args, char** argList, int16_t& prepros, string& filePref, int32_t& s, int32_t& k, int32_t& g, vector<string>& seqs, int32_t& t, string& qFile, string& c, uint32_t& m, SrchStrd& strd, bool& r, int16_t& X, uint16_t &nRes){
	int option_index = 0, a; 

	//Check wheather arguments are given for anything at all
	if(nb_args < MIN_PARAM_NB) return false;

	static struct option long_options[] = {
        {"graph-prefix",   required_argument,  0, 'i'},
        {"seed-length",    optional_argument,  0, 's'},
        {"kmer-length",    optional_argument,  0, 'k'},
        {"min-length",     optional_argument,  0, 'g'},
        {"fasta",          optional_argument,  0, 'f'},
        {"threads",        optional_argument,  0, 't'},
        {"query",          optional_argument,  0, 'q'},
        {"colors",         optional_argument,  0, 'c'},
        {"quorum",         optional_argument,  0, 'm'},
        {"X-dropoff",      optional_argument,  0, 'X'},
        {"max-results",    optional_argument,  0, 'n'},
        {"strand",         optional_argument,  0, 'd'},
        {"runtimes",       no_argument,        0, 'r'},
        {0,                0,                  0,  0 }
    };

    //Try to get the command
	if(!strcmp(argList[1], BUILD_COMMAND)){
		//Set preprocessing flag
		prepros = 0;
	} else if(!strcmp(argList[1], SEARCH_COMMAND)){
		//Set preprocessing flag
		prepros = -1;
	} else{
		cerr << "ERROR: Command unknown" << endl;
		return false;
	}

	//Testing
	// cout << "nb_args:" << nb_args << endl;
	// cout << "argList:" << endl;
	// for(int i = 0; i < nb_args; ++i){
	// 	cout << argList[i] << endl;
	// }
	// //a = getopt_long(nb_args, argList, OPTIONS, long_options, &option_index);
	// cout << "Do we get at least here?" << endl;

	//Parse all parameters given
	//while((a = getopt(nb_args, argList, OPTIONS)) != -1){
	while ((a = getopt_long(nb_args, argList, OPTIONS, long_options, &option_index)) != -1){
		//Testing
		// cout << "We get here\na:" << (char) a << endl;

		//Assign parameter values
		switch(a){
			case 'i':
				//Testing
				//cout << "optarg:" << optarg << endl;

				filePref = optarg;
				break;
			case 'f':
				//Testing
				/*cout << "Why do we not go here?" << endl;*/
				//cout << "optind:" << optind << endl;

				for(--optind; (optind < nb_args) && (*argList[optind] != '-'); ++optind){
					seqs.push_back(argList[optind]);
				}

				//Note that we have to build a graph
				++prepros;
				break;
			case 's':
				//Testing
				//cout << "Do we get here?" << endl;

				s = atoi(optarg);
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'g':
				g = atoi(optarg);
				break;
			case 't':
				t = atoi(optarg);
				break;
			case 'q':
				qFile = optarg;
				break;
			case 'c':
				c = optarg;
				break;
			case 'm':
				//A quorum has to be positive
				if(atoi(optarg) <= 0){
					cerr << "ERROR: Quorum parameter has to be a positive number" << endl;
					return false;
				} 
				
				m = atoi(optarg);
				break;
			case 'X':
				//Testing
				// cerr << atoi(optarg) << endl;

				//Check if X is applicable
				if(atoi(optarg) < 0 || atoi(optarg) > INT16_MAX){
					//Testing
					// if(atoi(optarg) < 0){
					// 	// cout << "6 Option 1" << endl;
					// } else{
					// 	cout << "6 Option 2" << endl;
					// }
					// if(atoi(optarg) > INT16_MAX){
					// 	cout << "7 Option 1" << endl;
					// } else{
					// 	// cout << "7 Option 2" << endl;
					// }

					cerr << "ERROR: X-dropoff value not applicable" << endl;
					return false;
				}

				X = atoi(optarg);
				break;
			case 'n':
				//Check if value is valid
				if(atoi(optarg) <= 0){
					//Testing
					// cout << "8 Option 1" << endl;

					cerr << "ERROR: Maximum number of alignments should be a positive number" << endl;
					return false;
				} else if(atoi(optarg) > UINT16_MAX){
					//Testing
					// cout << "9 Option 1" << endl;

					cerr << "Maximum number of alignments set to maximum" << endl;
					nRes = UINT16_MAX;
				} else{
					//Testing
					// cout << "9 Option 2" << endl;

					nRes = atoi(optarg);
				}

				//Testing
				// cout << "8 Option 2" << endl;

				break;
			case 'r':
				r = true;
				break;
			case 'd':
				//Testing
				// cout << optarg << endl;

				if(*optarg == PLUS_STRAND){
					strd = Plus;

					//Testing
					// cout << "10 Option 1" << endl;
					// exit(0);
				} else if(*optarg == MINUS_STRAND){
					//Testing
					// cout << "10 Option 2" << endl;

					strd = Minus;
				} else{
					//Testing
					// cout << "10 Option 3" << endl;
					
					cerr << "Unrecognized strand option" <<  endl << "Default is used" << endl;
				}
			default:
				//Testing
				//cout << "Default: optarg:" << optarg << endl;

				break;
		}
	}

	//Testing
	// cout << "a: " << a << endl;

	//Without a file prefix none command can be executed
	if(!strcmp(filePref.c_str(), "")) return false;

	//Testing
	//cout << "!strcmp(filePref.c_str(), ""):" << !strcmp(filePref.c_str(), "") << endl;
	
	//The search command additionally needs a query sequence
	if(prepros < 0 && !strcmp(qFile.c_str(), "")) return false;

	return true;
}

//This function reads in a file in which colors are stored the search will be based on
const list<pair<string, size_t>> loadSearchColors(const char* filename, uint32_t& nbCols){
	string line;
	pair<string, size_t> color("", -1);//Color id is set to -1 here to avoid false results if color from color search set is not present in the graph//TODO: Avoid this by mapping the colors directly within this function!
	list<pair<string, size_t>> colorlist;

	//Open the file
	ifstream fStr(filename);

	//Check if the file is open
	if(!fStr.is_open()){
		cerr << "ERROR: Color set file could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	//Reset counter
	nbCols = 0;
	//Read in lines of color set file
	while(getline(fStr, line)){
		//Store color name
		color.first = line;
		//Store each color
		colorlist.push_back(color);
		//Increment counter
		++nbCols;
	}

	//Close the file
	fStr.close();

	return colorlist;
}

//This function loads query sequences from a file and stores it in a list
void loadQueries(const string &filename, vector<string> &qList){
	string line;

	//Open the file
	ifstream fStr(filename);

	//Check if the file is open
	if(!fStr.is_open()){
		cerr << "ERROR: Query sequence file could not be opened" << endl;
		exit(EXIT_FAILURE);
	}

	//Read in queries (one per line)
	while(getline(fStr, line)) qList.push_back(line);
	//Close file
	fStr.close();
}