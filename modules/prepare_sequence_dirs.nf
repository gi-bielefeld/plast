process PREPARE_SEQUENCE_DIRS {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:d9fa884_2026-09-03"

    input:
    tuple path(refSeqFiles, arity: '0..*'),
          path(rawSeqFiles, arity: '0..*')

    output:
    tuple path('references'),
          path('reads'),
          emit: sequence_dirs

    script:
    """
    mkdir -p references
    mkdir -p reads

    if [ "\${#refSeqFiles[@]}" -gt 0 ]; then
        for input_file in "\${refSeqFiles[@]}"; do
            cp -a "\${input_file}" references/
        done
    fi

    if [ "\${#rawSeqFiles[@]}" -gt 0 ]; then
        for input_file in "\${rawSeqFiles[@]}"; do
            cp -a "\${input_file}" reads/
        done
    fi
    """
}
