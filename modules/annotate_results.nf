process ANNOTATE_RESULTS {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:7dc866f_2026-08-19"

    publishDir params.outdir, mode: 'copy'

    input:
    path(search_dir)

    tuple path(query_file), path(query_map_file), path(query_log_file)

    output:
    path 'plast_results'

    script:
    """
    mkdir -p plast_results

    python3 /plast/clowm/scripts/annotateResults.py \\
    "${search_dir}/results.plast" \\
    "${query_map_file}" \\
    "plast_results/results.annotated.plast"

    cp "${search_dir}/results.plast" plast_results/
    cp "${search_dir}/plast.stderr.txt" plast_results/
    cp "${search_dir}/search.command.txt" plast_results/
    """
}
