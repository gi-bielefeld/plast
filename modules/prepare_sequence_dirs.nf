process PREPARE_SEQUENCE_DIRS {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:d9fa884_2026-09-03"

    input:
    path refSeqFiles, arity: '0..*'
    path rawSeqFiles, arity: '0..*'

    output:
    tuple path('references'),
          path('reads'),
          emit: sequence_dirs

    script:
    """
    mkdir -p references
    mkdir -p reads

    if [ -z "${refSeqFiles}" ] && [ -z "${rawSeqFiles}" ]; then
        echo "No supported sequence files were found." >&2
        exit 1
    fi

    #if [ "\${#refSeqFiles[@]}" -gt 0 ]; then

    if [ -n "${refSeqFiles}" ]; then

        #for input_file in "\${refSeqFiles[@]}"; do

        for input_file in ${refSeqFiles}; do
            cp -a "\${input_file}" references/
        done
    else
        echo "No reference sequence files were provided"
    fi

    #if [ "\${#rawSeqFiles[@]}" -gt 0 ]; then

    if [ -n "${rawSeqFiles}" ]; then
        for input_file in ${rawSeqFiles}; do
            cp -a "\${input_file}" reads/
        done
    else
        echo "No read sequence files were provided"
    fi
    """
}
