process NORMALIZE_QUERY {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:7dc866f_2026-08-19"

    input:
      path query_fasta

    output:
      tuple path('*.q'), path('query.map.tsv'), path('query.normalization.log'), emit: normalized_query

    script:
    """
      python3 /plast/clowm/scripts/normalizeQuery.py "${query_fasta}" "${params.random_seed}"
    """
}
