process SEARCH_PLAST {
    label 'medium'
    container "ghcr.io/gi-bielefeld/plast:d9fa884_2026-09-03"

    //TODO: Is this necessary?
    // shell '/bin/bash', '-euo', 'pipefail'

    input:
    tuple path(graph_dir),
          path(graph_prefix_file),
          path(query_file),
          val(use_search_color_set),
          path(search_color_set_file)

          // path(query_map_file)

    output:
    path('search'), emit: search_results

    script:
    """
    mkdir -p search

    graph_prefix=\$(cat "${graph_prefix_file}")

    plast_args=(
        Search
        -i "${graph_dir}/\${graph_prefix}"
        -q "${query_file}"
        -w ${params.seed_length}
        -n ${params.max_results}
        -m ${params.mismatch}
        -M ${params.match}
        -d ${params.gap_score}
        -T ${params.evalue_threshold}
    )

    if [ "${params.strand}" != "both" ]; then
        plast_args+=( -o "${params.strand}" )
    fi

    if [ "${use_search_color_set}" = "true" ]; then
        plast_args+=( -s "${search_color_set_file}" )
    fi

    if [ "${params.quorum}" != "null" ]; then
        plast_args+=( -Q ${params.quorum} )
    fi

    if [ "${params.x_dropoff}" != "null" ]; then
        plast_args+=( -X ${params.x_dropoff} )
    fi

    if [ "${params.lambda_ungapped}" != "null" ]; then
        plast_args+=( -l ${params.lambda_ungapped} )
    fi

    if [ "${params.lambda_gapped}" != "null" ]; then
        plast_args+=( -L ${params.lambda_gapped} )
    fi

    if [ "${params.stat_c_ungapped}" != "null" ]; then
        plast_args+=( -c ${params.stat_c_ungapped} )
    fi

    if [ "${params.stat_c_gapped}" != "null" ]; then
        plast_args+=( -C ${params.stat_c_gapped} )
    fi

    if [ "${params.report_colors}" = "true" ]; then
        plast_args+=( -r )
    fi

    PLAST "\${plast_args[@]}" \\
        > search/results.plast \\
        2> search/plast.stderr.txt

    cp "${query_file}" search/
    printf 'PLAST' > search/search.command.txt

    for argument in "\${plast_args[@]}"; do
        printf ' %q' "\${argument}" >> search/search.command.txt
    done

    printf '\\n' >> search/search.command.txt
    """
}
