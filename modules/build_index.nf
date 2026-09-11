process BUILD_INDEX {
    label 'medium'
    container "ghcr.io/gi-bielefeld/plast:d9fa884_2026-09-03"

    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple path(graph_dir),
          path(graph_prefix_file)

    output:
    tuple path('indexed_graph'),
          path('indexed_graph/graph_prefix.txt'),
          emit: graph

    script:
    """
    graph_prefix=\$(cat "${graph_prefix_file}")

    mkdir -p indexed_graph
    cp -r "${graph_dir}/." indexed_graph/

    PLAST Build \\
        -i "indexed_graph/\${graph_prefix}" \\
        -w ${params.seed_length} \\
        > indexed_graph/build_index.stdout.txt \\
        2> indexed_graph/build_index.stderr.txt

    printf 'PLAST Build -i %q -w %q\\n' \\
        "indexed_graph/\${graph_prefix}" \\
        "${params.seed_length}" \\
        > indexed_graph/build_index.command.txt

    if [ ! -f "indexed_graph/\${graph_prefix}.idx" ]; then
        echo "PLAST Build did not create an index file." >&2
        exit 1
    fi

    cp "${graph_prefix_file}" indexed_graph/graph_prefix.txt
    """
}
