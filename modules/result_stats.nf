process RESULT_STATS {
    label 'tiny'

    container 'ghcr.io/gi-bielefeld/plast:d9fa884_2026-09-03'

    publishDir params.outdir, mode: 'copy'

    input:
    path search_dir

    // path stats_script

    output:
    path 'result_stats.txt',
         emit: result_stats

    script:
    """
    if [ ! -f "${search_dir}/results.plast" ]; then
        echo "PLAST result file was not found: ${search_dir}/results.plast" >&2
        exit 1
    fi

    python3 /plast/comparison/scripts/showPLASTresStats.py \
        "${search_dir}/results.plast" \
        > result_stats.txt
    """
}
