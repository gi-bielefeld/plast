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
        error "Parameter '--${parameter_name}' was not set."
    }

    if( is_s3_path(parameter_value) ) {
        return
    }

    def input_file = file(parameter_value)

    if( !input_file.exists() || input_file.isDirectory() ) {
        error "Parameter '--${parameter_name}' must point to an existing file: ${parameter_value}"
    }
}

def validate_input_directory(String parameter_name, String parameter_value) {
    if( !parameter_value ) {
        error "Parameter '--${parameter_name}' was not set."
    }

    if( is_s3_path(parameter_value) ) {
        return
    }

    def input_directory = file(parameter_value)

    if( !input_directory.exists() || !input_directory.isDirectory() ) {
        error "Parameter '--${parameter_name}' must point to an existing directory: ${parameter_value}"
    }
}

def validate_parameters() {
    if( !(params.mode in ['Build', 'Search', 'Build and Search']) ) {
        error "Parameter '--mode' must be 'Build', 'Search' or 'Build and Search'."
    }

    def kmer_length = params.kmer_length.toInteger()

    if( kmer_length < 9 ) {
        error "Parameter '--kmer_length' must be at least 9."
    }

    if( kmer_length > 63) {
        error "Parameter  '--kmer_length' must be at most 63."
    }

    def min_length = params.min_length.toInteger()

    if( min_length < 4 ) {
        error "Parameter '--min_length' must be must at least 4."
    }

    if( min_length > kmer_length - 2 ) {
        error "Parameter '--min_length' must not exceed --kmer_length - 2."
    }

    def seed_length = params.seed_length.toInteger()

    if( seed_length < 11 || seed_length > 15 ) {
        error "Parameter '--seed_length' must be between 11 and 15."
    }

    if( params.gap_score.toInteger() >= 0 ) {
        error "Parameter '--gap_score' must be negative."
    }

    if( params.quorum != null && params.quorum.toInteger() < 1 ) {
        error "Parameter '--quorum' must be at least 1."
    }

    if( params.strand !in ['+', '-', 'both'] ) {
        error "The parameter '--strand' must be '+', '-' or 'both'."
    }

    if( params.ref_seq_input_dir ) {
        validate_input_directory('ref_seq_input_dir', params.ref_seq_input_dir)
    }

    if( params.raw_seq_input_dir ) {
        validate_input_directory('raw_seq_input_dir', params.raw_seq_input_dir)
    }

    if( params.input_zip ) {
        validate_input_file('input_zip', params.input_zip)
    }

    if( params.input_graph_dir ) {
        if( !params.graph_prefix ) {
            error "Parameter '--graph_prefix' is required with --input_graph_dir."
        }

        if( !(
            params.input_graph_dir.startsWith('s3://') ||
            params.input_graph_dir.startsWith('/') ||
            params.input_graph_dir.startsWith('./') ||
            params.input_graph_dir.startsWith('../')
        ) ) {
            error "Parameter '--input_graph_dir' must be a local path or an s3:// URI."
        }
    }

    if( params.search_color_set ) {
        validate_input_file('search_color_set', params.search_color_set)
    }

    def query_required = params.mode in ['Search', 'Build and Search']

    if( query_required && !params.query_fasta ) {
        error "A query file is needed to perform a PLAST search"
    }
    
    if ( query_required ) {
        validate_input_file('query_fasta', params.query_fasta)
    }
}

def sequence_pattern(String directory) {
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

    if( params.input_graph_dir ) {
        has_existing_graph = true
    }

    if( params.mode == 'Search' ) {
        if( !has_existing_graph ) {
            error(
                "Mode 'Search' requires an existing graph with both " +
                ".gfa/.gfa.gz and .color.bfg files."
            )
        }
    }

    if( params.mode == 'Build' && !has_existing_graph && !has_sequence_directories ) {
        error(
            "Mode 'Build' requires either an existing graph or input sequences."
        )
    }

    if( params.mode == 'Build and Search' &&
        !has_existing_graph &&
        !has_sequence_directories ) {
        error(
            "Mode 'Build and Search' requires either an existing graph or input sequences."
        )
    }

    log.info "Starting PLAST workflow"
    log.info "Mode: ${params.mode}"

    if( params.query_fasta ) {
        log.info "Query file: ${params.query_fasta}"
    }
    else {
        log.info "No query file was provided. The workflow will end after graph/index construction."
    }

    log.info "Output directory: ${params.outdir}"

    //If sequences for graph construction are given, a new graph will always be built in mode "Build and Search"
    def build_new_graph_from_sequences =
        (params.mode == 'Build and Search' && has_sequence_directories) ||
        (params.mode == 'Build' && !has_existing_graph && has_sequence_directories)

    def build_graph_prefix = params.graph_prefix ?: 'plast_graph'

    if( !build_new_graph_from_sequences && !params.graph_prefix ) {
        error(
            "Parameter '--graph_prefix' is required if no new graph is built."
        )
    }

    if( has_existing_graph &&
        has_sequence_directories &&
        !build_new_graph_from_sequences ) {
        log.info(
            "Both an existing graph and sequence directories were provided. " +
            "The existing graph will be used; sequence directories are ignored."
        )
    }

    if( params.input_graph_dir && params.graph_prefix ){
        graph_pattern = "${params.input_graph_dir}/${params.graph_prefix}.*"

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
        if( params.ref_seq_input_dir ) {
            ref_files_ch = channel
                .fromPath(sequence_pattern(params.ref_seq_input_dir), checkIfExists: true)
                .filter { sequence_file(it) }
                .collect()
        } else{
            ref_files_ch = channel.of([])
        }

        if( params.raw_seq_input_dir ){
            raw_files_ch = channel
                .fromPath(sequence_pattern(params.raw_seq_input_dir), checkIfExists: true)
                .filter { sequence_file(it) }
                .collect()
        } else{
            raw_files_ch = channel.of([])
        }

        // sequence_file_lists_ch = ref_files_ch
        //     .collect()
        //     .combine(raw_files_ch.collect())
        //     .map { reference_files, read_files ->
        //         tuple(reference_files, read_files)
        //     }

        sequence_dirs_ch = PREPARE_SEQUENCE_DIRS(ref_files_ch, raw_files_ch)
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
        log.info "Building a new PLAST graph from sequence input."
        graph_ch = BUILD_GRAPH(sequence_dirs_ch)
    }
    else if( params.mode == 'Build' || params.mode == 'Build and Search' ) {
        log.info "Building a new PLAST index for the existing graph."
        graph_ch = BUILD_INDEX(existing_graph_ch)
    }
    else {
        graph_ch = existing_graph_ch
    }

    //If mode is "Build", we are done
    if( params.mode == 'Build' ) {
        return
    }

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
        search_color_set_ch = channel
            .fromPath(
                params.search_color_set,
                checkIfExists: true
            )
            .map { search_color_set_file ->
                tuple(true, search_color_set_file)
            }
    } else{
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

    RESULT_STATS(SEARCH_PLAST.out.search_results)

    ANNOTATE_RESULTS(SEARCH_PLAST.out.search_results, normalized_query_ch)
}
