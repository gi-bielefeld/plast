include {
    NORMALIZE_QUERY
} from './modules/normalize_query'

include {
    PREPARE_SEQUENCE_DIRS
} from './modules/prepare_sequence_dirs'

include {
    PREPARE_EXISTING_GRAPH
} from './modules/prepare_existing_graph'

include {
    PREPARE_ZIP_GRAPH
} from './modules/prepare_zip_graph'

include {
    BUILD_GRAPH
} from './modules/build_graph'

include {
    BUILD_INDEX
} from './modules/build_index'

include {
    SEARCH_PLAST
} from './modules/search_plast'

include {
    RESULT_STATS
} from './modules/result_stats'

include {
    ANNOTATE_RESULTS
} from './modules/annotate_results'

def is_s3_path(String value) {
    return value != null && value.startsWith('s3://')
}

def validate_input_file(String parameter_name, String parameter_value) {
    if( !parameter_value ) {
        //Testing
        log.info "Input file parameter was not set"

        error "Parameter '--${parameter_name}' was not set."
    }

    if( is_s3_path(parameter_value) ) {
        //Testing
        log.info "Tests 3 and 7: A given file parameter value is an s3 path"

        return
    }

    //Testing
    log.info "Tests 2, 4, 5, 6, 7, and 8: A given parameter value is a conventional file path"

    def input_file = file(parameter_value)

    //Testing
    if( !input_file.exists() ){
        log.info "Specified input file does not exist"
    }
    if( input_file.isDirectory() ){
        log.info "Specified input file is a directory"
    }

    if( !input_file.exists() || input_file.isDirectory() ) {
        error "Parameter '--${parameter_name}' must point to an existing file: ${parameter_value}"
    }

    //Testing
    log.info "Tests 2, 4, 5, 6, 7, and 8: Specified input file exists"
}

def validate_input_directory(String parameter_name, String parameter_value) {
    if( !parameter_value ) {
        //Testing
        log.info "Input directory parameter was not set"

        error "Parameter '--${parameter_name}' was not set."
    }

    if( is_s3_path(parameter_value) ) {
        //Testing
        log.info "A given directory parameter value is an s3 path"

        return
    }

    //Testing
    log.info "A given parameter value is a conventional directory path"

    def input_directory = file(parameter_value)

    //Testing
    if( !input_directory.exists() ){
        log.info "Specified input directory does not exist"
    }
    if( !input_directory.isDirectory() ){
        log.info "Specified input directory is not a directory"
    }

    if( !input_directory.exists() || !input_directory.isDirectory() ) {
        error "Parameter '--${parameter_name}' must point to an existing directory: ${parameter_value}"
    }

    //Testing
    log.info "Test 1: Specified input directory exists"
}

// def graph_files(String graph_dir, String prefix) {
//     def directory = file(graph_dir)

//     //Testing
//     log.info "${directory}/${prefix}.gfa"

//     return [
//         gfa: file("${directory}/${prefix}.gfa"),
//         gfa_gz: file("${directory}/${prefix}.gfa.gz"),
//         colors: file("${directory}/${prefix}.color.bfg"),
//         index: file("${directory}/${prefix}.idx")
//     ]
//}

def validate_parameters() {
    if( !(params.mode in ['Build', 'Search', 'Build and Search']) ) {
        //Testing
        log.info "Given mode is not one of the expected"

        error "Parameter '--mode' must be 'Build', 'Search' or 'Build and Search'."
    }

    //Testing
    log.info "Tests 1, 2, 3, 4, 5, 6, 7, and 8: Mode parameter is one of the expected"

    // def modes = [
    //     params.graph_gfa ? 'gfa' : null,
    //     params.reference_dir && params.reads_dir ? 'directories' : null,
    //     params.graph_zip ? 'zip' : null
    // ].findAll { it != null }

    // if( modes.size() != 1 ) {
    //     error(
    //         "Exactly on graph input option must be chosen: " +
    //         "--graph_gfa or --reference_dir together with --reads_dir or --graph_zip."
    //     )
    // }

    def kmer_length = params.kmer_length.toInteger()

    if( kmer_length < 9 ) {
        //Testing
        log.info "Given k is <9"

        error "Parameter '--kmer_length' must be at least 9."
    } else{
        //Testing
        log.info "Tests 1, 2, 3, 4, 5, 6, 7, and 8: Given k is >= 9"
    }

    if( kmer_length > 63) {
        //Testing
        log.info "Given k >63"

        error "Parameter  '--kmer_length' must be at most 63."
    } else{
        //Testing
        log.info "Tests 1, 2, 3, 4, 5, 6, 7, and 8: Given k <=63"
    }

    def min_length = params.min_length.toInteger()

    if( min_length < 4 ) {
        //Testing
        log.info "Given minimizer length is <4"

        error "Parameter '--min_length' must be must at least 4."
    } else{
        //Testing
        log.info "Tests 1, 2, 3, 4, 5, 6, 7, and 8: Given minimizer length is >=4"
    }

    if( min_length > kmer_length - 2 ) {
        //Testing
        log.info "Given minimizer length >k-2"

        error "Parameter '--min_length' must not exceed --kmer_length - 2."
    } else{
        //Testing
        log.info "Tests 1, 2, 3, 4, 5, 6, 7, and 8: Given minimizer length is <=k-2"
    }

    def seed_length = params.seed_length.toInteger()

    if( seed_length < 11 || seed_length > 15 ) {
        //Testing
        if( seed_length < 11 ){
            log.info "Given seed length is <11"
        }
        if( seed_length > 15){
            log.info "Given seed length is >15"
        }

        error "Parameter '--seed_length' must be between 11 and 15."
    } else{
        //Testing
        log.info "Test 1, 2, 3, 4, 5, 6, 7, and 8: Given seed length is between 11 and 15"
    }

    if( params.gap_score.toInteger() >= 0 ) {
        //Testing
        log.info "Given gap score is >=0"

        error "Parameter '--gap_score' must be negative."
    } else{
        //Testing
        log.info "Tests 1, 2, 3, 4, 5, 6, 7, and 8: Given gap score is <0"
    }

    if( params.quorum != null && params.quorum.toInteger() < 1 ) {
        //Testing
        log.info "Given quorum is <1"

        error "Parameter '--quorum' must be at least 1."
    } else{
        //Testing
        if( params.quorum != null ){
            log.info "Test 7: Given quorum is >=1"
        }
    }

    if( params.strand !in ['+', '-', 'both'] ) {
        //Testing
        log.info "Given value for parameter strand is unknown"

        error "The parameter '--strand' must be '+', '-' or 'both'."
    } else{
        //Testing
        log.info "Tests 1, 2, 3, 4, 5, 6, 7, and 8: Given value for parameter strand known"
    }

    if( params.ref_seq_input_dir ) {
        //Testing
        log.info "Parameter ref_seq_input_dir is set"

        validate_input_directory('ref_seq_input_dir', params.ref_seq_input_dir)
    }

    if( params.raw_seq_input_dir ) {
        //Testing
        log.info "Parameter raw_seq_input_dir is set"

        validate_input_directory('raw_seq_input_dir', params.raw_seq_input_dir)
    }

    if( params.input_zip ) {
        //Testing
        log.info "Test 8: Parameter input_zip is set"

        validate_input_file('input_zip', params.input_zip)
    }

    if( params.input_graph_dir ) {
        //Testing
        // log.info "Tests 1 and 2: Parameter input_graph_dir is set"

        // validate_input_directory('input_graph_dir', params.input_graph_dir)

        if( !params.graph_prefix ) {
            //Testing
            log.info "Parameter input_graph_dir is set, but no graph prefix is given"

            error "Parameter '--graph_prefix' is required with --input_graph_dir."
        } else{
            //Testing
            log.info "Tests 1, 2, 3, 4, 5, 6, and 7: Parameters input_graph_dir and graph_prefix are set"
        }

        if( !(
            params.input_graph_dir.startsWith('s3://') ||
            params.input_graph_dir.startsWith('/') ||
            params.input_graph_dir.startsWith('./') ||
            params.input_graph_dir.startsWith('../')
        ) ) {
            //Testing
            log.info "Parameter input_graph_dir is not local path or an s3 URI"

            error "Parameter '--input_graph_dir' must be a local path or an s3:// URI."
        }
    }

    // if( params.graph_prefix && !params.input_graph_dir ) {
    //     error "Parameter '--input_graph_dir' is required with --graph_prefix."
    // }

    // def search_color_set_ch = channel.empty()

    if( params.search_color_set ) {
        //Testing
        log.info "Tests 3 and 7: Parameter search_color_set is given"

        validate_input_file('search_color_set', params.search_color_set)

    //     search_color_set_ch = channel.fromPath(
    //         params.search_color_set,
    //         checkIfExists: true
    //     )

    }

    def query_required = params.mode in ['Search', 'Build and Search']

    if( query_required && !params.query_fasta ) {
        //Testing
        log.info "Workflow mode involves a search, but parameter query_fasta is not set"

        error "A query file is needed to perform a PLAST search"
    }
    
    //Testing    
    if ( !query_required && params.query_fasta ) {
        log.info "Workflow mode does not involve a search, but parameter query_fasta is set"
    }
    if( !query_required && !params.query_fasta ){
        log.info "Test 1: Workflow mode does not involve a search and parameter query_fasta is not set"
    }

    if ( query_required ) {
        //Testing
        log.info "Tests 2, 3, 4, 5, 6, 7, and 8: Workflow mode involves a search and parameter query_fasta is set"

        validate_input_file('query_fasta', params.query_fasta)
    }
}

def sequence_pattern(String directory) {
    //Testing
    if( directory.endsWith('/') ){
        log.info "Given input directory ends with /"
    } else{
        log.info "Given input directory does not end with /"
    }

    return directory.endsWith('/') ? "${directory}*" : "${directory}/*"
}

def sequence_file(path) {
    def name = path.name.toLowerCase()

    return name.endsWith('.fa') ||
           name.endsWith('.fasta') ||
           name.endsWith('.fna') ||
           name.endsWith('.fa.gz') ||
           name.endsWith('.fasta.gz') ||
           name.endsWith('.fna.gz') ||
           name.endsWith('.fof')
}

workflow {
    validate_parameters()

    def has_sequence_directories =
        params.ref_seq_input_dir ||
        params.raw_seq_input_dir ||
        params.input_zip

    def has_existing_graph = false

    // def has_graph_index = false
    // def has_partial_graph = false

    if( params.input_graph_dir ) {
        //Testing
        log.info "Tests 1, 2, 3, 4, 5, 6, and 7: Parameter input_graph_dir is given"
    //     log.info "params.input_graph_dir: ${params.input_graph_dir}"

    //     def graph = graph_files(
    //         params.input_graph_dir,
    //         params.graph_prefix
    //     )

    //     def has_gfa = graph.gfa.exists() || graph.gfa_gz.exists()
    //     def has_colors = graph.colors.exists()

        has_existing_graph = true

    //     //Testing
    //     if (has_gfa) {
    //         log.info "GFA file has been found"
    //     }
    //     if(has_colors) {
    //         log.info "Color file has been found"
    //     }

    //     has_graph_index = graph.index.exists()
    //     has_partial_graph = has_gfa != has_colors

    }

    // if( has_partial_graph && !has_sequence_directories ) {
    //     //Testing
    //     log.info "Not all required graph files are present and no sequence directories are given"

    //     error(
    //         "A graph requires both a .gfa/.gfa.gz file and a .color.bfg file. " +
    //         "Only one of these files is present and no sequence input was provided."
    //     )
    // } else{
    //     //Testing
    //     if ( has_partial_graph ) {
    //         log.info "Not all required graph files are present and sequence directories are given"
    //     } else{
    //         if ( has_sequence_directories ){
    //             log.info "A graph is not partially present and sequence directories are given"
    //         } else{
    //             log.info "Test 1: A graph is not partially present and sequence directories are not given"
    //         }
    //     }
    // }

    if( params.mode == 'Search' ) {
        //Testing
        log.info "Tests 2, 3, 4, 5, 6, and 7: Workflow mode is 'Search'"

        if( !has_existing_graph ) {
            //Testing
            log.info "Workflow mode is search and there is no existing graph"

            error(
                "Mode 'Search' requires an existing graph with both " +
                ".gfa/.gfa.gz and .color.bfg files."
            )
        }

        //Testing
        log.info "Tests 2, 3, 4, 5, 6, and 7: Workflow mode is Search and there is an existing graph"

        // if( !has_graph_index ) {
        //     //Testing
        //     log.info "Workflow mode is Search and there is no index file"

        //     error(
        //         "Mode 'Search' requires an existing .idx file."
        //     )
        // }

        // //Testing
        // log.info "Workflow mode is Search and all required graph files exist"
    }

    if( params.mode == 'Build' && !has_existing_graph && !has_sequence_directories ) {
        //Testing
        log.info "Workflow mode is Build and there is neither an existing graph nor sequence directories"

        error(
            "Mode 'Build' requires either an existing graph or input sequences."
        )
    } else{
        //Testing
        if( params.mode == 'Build' && !has_existing_graph ){
            log.info "Workflow mode is Build and there is no existing graph, but sequence directories"
        }
        if( params.mode == 'Build' && !has_sequence_directories ){
            log.info "Test 1: Workflow mode is Build and there are no sequence directories, but an existing graph"
        }
        if( !has_existing_graph && !has_sequence_directories ) {
            log.info "There is an existing graph and sequence directories, but workflow mode is not Search"
        }
    }

    //Testing
    if( params.mode == 'Build' && has_existing_graph && has_sequence_directories ){
        log.info "Workflow mode is Build, there is an existing graph and there are sequence directories"
    }
    if( params.mode != 'Build' && has_existing_graph && !has_sequence_directories ){
        log.info(
            "Tests 2, 3, 4, 5, 6, and 7: Workflow mode is not Build, there is an existing graph, but no sequence " +
            "directories"
        )
    }
    if( params.mode != 'Build' && !has_existing_graph && has_sequence_directories){
        log.info "Test 8: Workflow mode is not Build, there is no graph and there are sequence directories"
    }
    if( params.mode != 'Build' && !has_existing_graph && !has_sequence_directories){
        log.info "Workflow mode is not Build, there is no graph and no sequence directories"
    }

    if( params.mode == 'Build and Search' &&
        !has_existing_graph &&
        !has_sequence_directories ) {
        //Testing
        log.info "Workflow mode is Build and Search, but there is not graph and no sequence directories"

        error(
            "Mode 'Build and Search' requires either an existing graph or input sequences."
        )
    }

    //Testing
    if ( params.mode == 'Build and Search' && !has_existing_graph && has_sequence_directories ){
        log.info "Test 8: Workflow mode is Build and Search, no graph, but sequence directories"
    }
    if ( params.mode == 'Build and Search' && has_existing_graph && !has_sequence_directories ){
        log.info "Workflow is Build and Search, graph exists, but no sequence directories"
    }
    if ( params.mode == 'Build and Search' && has_existing_graph && has_sequence_directories ){
        log.info "Workflow is Build and Search graph exists and sequence directories exist"
    }
    if ( params.mode != 'Build and Search' && !has_existing_graph && has_sequence_directories ){
        log.info "Workflow mode is not Build and Search, no graph, but sequence directories"
    }
    if ( params.mode != 'Build and Search' && has_existing_graph && !has_sequence_directories ){
        log.info(
            "Tests 1, 2, 3, 4, 5, 6, and 7: Workflow mode is not Build and Search, graph exists, but no sequence " +
            "directories"
        )
    }
    if ( params.mode != 'Build and Search' && has_existing_graph && has_sequence_directories ){
        log.info "Workflow mode is Build and Search, graph exists and sequence directories exist"
    }
    if ( params.mode != 'Build and Search' && !has_existing_graph && !has_sequence_directories ){
        log.info "Workflow mode is not Build and Search, no graph and no sequence directories"
    }

    log.info "Starting PLAST workflow"
    log.info "Mode: ${params.mode}"

    if( params.query_fasta ) {
        //Testing
        log.info "Tests 2, 3, 4, 5, 6, 7, and 8: Query file is given"

        log.info "Query file: ${params.query_fasta}"
    }
    else {
        //Testing
        log.info "Test 1: No query file is given"

        log.info "No query file was provided. The workflow will end after graph/index construction."
    }

    log.info "Output directory: ${params.outdir}"

    //If sequences for graph construction are given, a new graph will always be built in mode "Build and Search"
    def build_new_graph_from_sequences =
        (params.mode == 'Build and Search' && has_sequence_directories) ||
        (params.mode == 'Build' && !has_existing_graph && has_sequence_directories)

    def build_graph_prefix = params.graph_prefix ?: 'plast_graph'

    if( !build_new_graph_from_sequences && !params.graph_prefix ) {
        //Testing
        log.info "No new graph shall be built and no graph prefix is given"

        error(
            "Parameter '--graph_prefix' is required if no new graph is built."
        )
    }

    //Testing
    if ( !build_new_graph_from_sequences && params.graph_prefix ) {
        log.info "Tests 1, 2, 3, 4, 5, 6, and 7: No new graph shall be built and graph_prefix is given"
    }
    if ( build_new_graph_from_sequences && !params.graph_prefix ){
        log.info "New graph shall be built and no graph prefix is given"
    }
    if ( build_new_graph_from_sequences && params.graph_prefix ){
        log.info "Test 8: New graph shall be built and graph_prefix is given"
    }

    // if( has_partial_graph && has_sequence_directories ) {
    //     //Testing
    //     log.info "Graph is only partially present and sequence directories are present"

    //     log.info(
    //         "Incomplete graph input detected. A new graph will be built from input sequences."
    //     )
    // }

    if( has_existing_graph &&
        has_sequence_directories &&
        !build_new_graph_from_sequences ) {
        //Testing
        log.info "A graph and sequences are provided"

        log.info(
            "Both an existing graph and sequence directories were provided. " +
            "The existing graph will be used; sequence directories are ignored."
        )
    }

    // def graph_input_description =
    // params.graph_gfa ? 'existing graph' :
    // params.graph_zip ? 'zip archive' :
    // 'sequence directory'

    // log.info "Graph input type: ${graph_input_description}"
    // log.info "Query file: ${params.query_fasta}"
    // log.info "Output directory: ${params.outdir}"

    //Prepare existing graph

    // if( has_existing_graph ) {

    if( params.input_graph_dir && params.graph_prefix ){
        //Testing
        // log.info "Test 1: An existing graph is provided"
        // log.info "Parameters input_graph_dir and graph_prefix are both set"

        graph_pattern = "${params.input_graph_dir}/${params.graph_prefix}.*"

        // existing_graph_input_ch = channel.of(
        //     tuple(
        //         params.graph_prefix,
        //         file(params.input_graph_dir)
        //     )
        // )

        existing_graph_files_ch = channel
            .fromPath(graph_pattern, checkIfExists: true)
            .collect()
            .map { graph_files ->
                tuple(params.graph_prefix, graph_files)
            }

        existing_graph_ch = PREPARE_EXISTING_GRAPH(existing_graph_files_ch)
    }

    //Prepare graph from zip archive
    if( params.input_zip ) {
        //Testing
        log.info "Test 8: A ZIP archive is provided"

        zip_input_ch = channel
            .fromPath(params.input_zip, checkIfExists: true)

        sequence_dirs_ch = PREPARE_ZIP_GRAPH(zip_input_ch)
            .map{ reference_dir, reads_dir ->
                tuple(
                    build_graph_prefix,
                    reference_dir,
                    reads_dir
                )
            }
    }

    //Prepare graph from reference and/or read directories
    else if( params.ref_seq_input_dir || params.raw_seq_input_dir ) {
        
        // def sequence_dirs = []

        if( params.ref_seq_input_dir ) {
            //Testing
            log.info "Reference sequences are provided"

            // sequence_dirs << file(params.reference_dir)

            ref_files_ch = channel
                .fromPath(sequence_pattern(params.ref_seq_input_dir), checkIfExists: true)
                .filter { sequence_file(it) }
        } else{
            ref_files_ch = channel.empty()
        }

        if( params.raw_seq_input_dir ){
            //Testing
            log.info "Read sequences are provided"

            // sequence_dirs << file(params.reads_dir) 

            raw_files_ch = channel
                .fromPath(sequence_pattern(params.raw_seq_input_dir), checkIfExists: true)
                .filter { sequence_file(it) }
        } else{
            raw_files_ch = channel.empty()
        }

        // direct_sequence_dirs_ch = channel.of(
        //     tuple(
        //         params.reference_dir ?
        //             file(params.reference_dir).name :
        //             '',
        //         params.reads_dir ?
        //             file(params.reads_dir).name :
        //             '',
        //         sequence_dirs
        //     )
        // )

        sequence_file_lists_ch = ref_files_ch
            .collect()
            .combine(raw_files_ch.collect())
            .map { reference_files, read_files ->
                tuple(reference_files, read_files)
            }

        sequence_dirs_ch = PREPARE_SEQUENCE_DIRS(sequence_file_lists_ch)
            .map{ reference_dir, reads_dir ->
                tuple(
                    build_graph_prefix,
                    reference_dir,
                    reads_dir
                )
            }
    }

    //Build a new graph and index if applicable
    if( build_new_graph_from_sequences ) {
        //Testing
        log.info "Test 8: A new graph is built"

        log.info "Building a new PLAST graph from sequence input."
        graph_ch = BUILD_GRAPH(sequence_dirs_ch)
    }
    else if( params.mode == 'Build' || params.mode == 'Build and Search' ) {
        //Testing
        log.info "Test 1: A new index is built"

        log.info "Building a new PLAST index for the existing graph."
        graph_ch = BUILD_INDEX(existing_graph_ch)
    }
    // else if( params.mode == 'Search' && params.seed_length > 11 ) {
    //     log.info(
    //         "seed_length > 11. A new PLAST index will be built before searching."
    //     )
    //     graph_ch = BUILD_INDEX(existing_graph_ch)
    // }
    else {
        //Testing
        log.info "Tests 2, 3, 4, 5, 6, and 7: An existing graph and index are used"

        graph_ch = existing_graph_ch
    }

    //If mode is "Build", we are done
    if( params.mode == 'Build' ) {
        //Testing
        log.info "Test 1: Workflow mode is Build"

        return
    }

    //Testing
    log.info "Tests 2, 3, 4, 5, 6, 7, and 8: Workflow mode is not Build"

    /*
     * We first need to transform the FASTA query file into the PLAST compatible format
     */
    query_ch = channel
        .fromPath(params.query_fasta, checkIfExists: true)

    normalized_query_ch = NORMALIZE_QUERY(query_ch)
    query_for_search_ch = normalized_query_ch.map {
        query_file,
        query_map_file,
        query_log_file ->

        query_file
    }

    if( params.search_color_set ){
        //Testing
        log.info "Tests 3 and 7: A search color set is provided"

        search_color_set_ch = channel
            .fromPath(
                params.search_color_set,
                checkIfExists: true
            )
            .map { search_color_set_file ->
                tuple(true, search_color_set_file)
            }
    } else{
        //Testing
        log.info "Tests 2, 4, 5, 6, and 8: No search color set is provided"

        search_color_set_ch = channel
            .fromPath(
                "${projectDir}/clowm/empty_search_color_set.txt",
                checkIfExists: true
            )
            .map { empty_search_color_set_file ->
                tuple(false, empty_search_color_set_file)
            }
    }

    search_ch = graph_ch
        .combine(query_for_search_ch)
        .combine(search_color_set_ch)

    SEARCH_PLAST(search_ch)

    // SEARCH_PLAST(graph_ch, query_for_search_ch)

    RESULT_STATS(SEARCH_PLAST.out.search_results)

    // annotated_result_ch = ANNOTATE_RESULTS(result_ch)

    ANNOTATE_RESULTS(SEARCH_PLAST.out.search_results, normalized_query_ch)
}
