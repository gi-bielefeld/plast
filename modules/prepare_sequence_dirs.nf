process PREPARE_SEQUENCE_DIRS {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:7dc866f_2026-08-19"

    input:
    tuple path(refSeqFiles, arity: '0..*'),
          path(rawSeqFiles, arity: '0..*')

    // tuple val(reference_dir_name),
    //       val(reads_dir_name),
    //       path(sequence_dirs)

    output:
    tuple path('references'),
          path('reads'),
          emit: sequence_dirs

    script:
    """
    mkdir -p references
    mkdir -p reads

    if [ "\${#refSeqFiles[@]}" -gt 0 ]; then
        #Testing
        echo Copying reference sequence files

        for input_file in "\${refSeqFiles[@]}"; do
            cp -a "\${input_file}" references/
        done
    fi

    if [ "\${#rawSeqFiles[@]}" -gt 0 ]; then
        #Testing
        echo Copying read sequence files

        for input_file in "\${rawSeqFiles[@]}"; do
            cp -a "\${input_file}" reads/
        done
    fi
    """
}
