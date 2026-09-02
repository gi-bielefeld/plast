process PREPARE_ZIP_GRAPH {
    label 'tiny'
    container "ghcr.io/gi-bielefeld/plast:7dc866f_2026-08-19"

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
        echo "The ZIP archive contains neither a 'references' nor a 'reads' directory." >&2
        exit 1
    fi

    if [ -n "\${reference_dir}" ]; then
        cp -a "\${reference_dir}/." references/
    fi

    if [ -n "\${reads_dir}" ]; then
        cp -a "\${reads_dir}/." reads/
    fi
    """
}
