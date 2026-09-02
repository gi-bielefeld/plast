process PREPARE_ZIP_GRAPH {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:7dc866f_2026-08-19"

    //TODO: Is this necessary?
    // shell '/bin/bash', '-euo', 'pipefail'

    input:
    path input_zip

    output:
    tuple path('references'), path('reads'), emit: sequence_dirs

    script:
    """
    mkdir -p extracted
    mkdir -p references
    mkdir -p reads

    unzip -q "${input_zip}" -d extracted

    reference_dir=\$(find extracted -type d -name references -print -quit)
    reads_dir=\$(find extracted -type d -name reads -print -quit)

    if [ -z "\${reference_dir}" ] && [ -z "\${reads_dir}" ]; then
        #Testing
        echo No sequence folder is found in ZIP archive

        echo "The ZIP archive contains neither a 'references' nor a 'reads' directory." >&2
        exit 1
    fi

    if [ -n "\${reference_dir}" ]; then
        #Testing
        echo Test 8: A reference sequence folder is found in ZIP archive

        cp -a "\${reference_dir}/." references/
    fi

    if [ -n "\${reads_dir}" ]; then
        #Testing
        echo A read folder is found in ZIP archive

        cp -a "\${reads_dir}/." reads/
    fi
    """
}
