process PREPARE_EXISTING_GRAPH {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:7dc866f_2026-08-19"

    //TODO: Is this necessary?
    // shell '/bin/bash', '-euo', 'pipefail'

    input:
    tuple val(graph_prefix),
          path(graph_files)

          // path(input_graph_dir)

    output:
    tuple path('graph'),
          path('graph/graph_prefix.txt'),
          emit: graph

    script:
    // def graph_base = graph_gfa.name
    //     .replaceFirst(/\.gz$/, '')
    //     .replaceFirst(/\.gfa$/, '')

    """
    mkdir -p graph

    for graph_file in ${graph_files}; do
        cp "\${graph_file}" graph/
    done

    if [ -f "graph/${graph_prefix}.gfa" ]; then
        #Testing
        echo Tests 1, 2, 3, 4, 5, 6, and 7: Uncompressed GFA file is found
    elif [ -f "graph/${graph_prefix}.gfa.gz" ]; then
        #Testing
        echo Compressed GFA file is found
    else
        #Testing
        echo No GFA file is found for the given graph prefix

        echo "No GFA file found for graph prefix '${graph_prefix}'." >&2
        exit 1
    fi

    if [ ! -f "graph/${graph_prefix}.color.bfg" ]; then
        #Testing
        echo No graph color file can be found

        echo "No .color.bfg file found for graph prefix '${graph_prefix}'." >&2
        exit 1
    fi

    #Testing
    echo Tests 1, 2, 3, 4, 5, 6, and 7: A graph color file is found

    if [ "${params.mode}" = "Search" ] &&
       [ ! -f "graph/${graph_prefix}.idx" ]; then
        #Testing
        echo Workflow mode is Search and no index file exists

        echo "Mode 'Search' requires '${graph_prefix}.idx'." >&2
        exit 1
    fi

    #Testing
    if [ "${params.mode}" = "Search" ]; then
        echo Tests 2, 3, 4, 5, 6, and 7: Workflow mode is Search and index file exists
    fi

    printf '%s\\n' "${graph_prefix}" > graph/graph_prefix.txt
    """
}
