process PREPARE_EXISTING_GRAPH {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:7dc866f_2026-08-19"

    input:
    tuple val(graph_prefix),
          path(graph_files)

    output:
    tuple path('graph'),
          path('graph/graph_prefix.txt'),
          emit: graph

    script:
    """
    mkdir -p graph

    for graph_file in ${graph_files}; do
        cp "\${graph_file}" graph/
    done

    if [ ! -f "graph/${graph_prefix}.gfa" ] && [ ! -f "graph/${graph_prefix}.gfa.gz" ]; then
        echo "No GFA file found for graph prefix '${graph_prefix}'." >&2
        exit 1
    fi

    if [ ! -f "graph/${graph_prefix}.color.bfg" ]; then
        echo "No .color.bfg file found for graph prefix '${graph_prefix}'." >&2
        exit 1
    fi

    if [ "${params.mode}" = "Search" ] &&
       [ ! -f "graph/${graph_prefix}.idx" ]; then
        echo "Mode 'Search' requires '${graph_prefix}.idx'." >&2
        exit 1
    fi

    printf '%s\\n' "${graph_prefix}" > graph/graph_prefix.txt
    """
}
