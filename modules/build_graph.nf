process BUILD_GRAPH {
    label 'medium'
    container "ghcr.io/gi-bielefeld/plast:d9fa884_2026-09-03"

    publishDir params.outdir, mode: 'copy'

    //TODO: Is this necessary?
    // shell '/bin/bash', '-euo', 'pipefail'

    input:
    tuple val(graph_prefix), 
          path(reference_dir), 
          path(reads_dir)

    output:
    tuple path('graph'),
          path('graph/graph_prefix.txt'),
          emit: graph

    script:
    """
    mkdir -p graph

    resolve_inputs() {
        local input_dir="\$1"
        local fof_file=\$(find "\${input_dir}" -maxdepth 1 -type f -name '*.fof' -print -quit)
        local file_name
        local resolved_file

        if [ -n "\${fof_file}" ]; then
            while IFS= read -r file_name || [ -n "\${file_name}" ]; do
                file_name="\$(echo "\${file_name}" | sed 's/^[[:space:]]*//;s/[[:space:]]*\$//')"

                if [ -z "\${file_name}" ]; then
                    continue
                fi

                resolved_file="\${input_dir}/\${file_name}"

                if [ ! -f "\${resolved_file}" ]; then
                    echo "The file listed in \${fof_file} does not exist: \${resolved_file}" >&2
                    exit 1
                fi

                printf '%s\\n' "\${resolved_file}"
            done < "\${fof_file}"
        else
            find -L "\${input_dir}" -maxdepth 1 -type f \\
                \\( -iname "*.fa" -o -iname "*.fna" -o -iname "*.fasta" \\
                -o -iname "*.fq" -o -iname "*.fastq" \\
                -o -iname "*.fa.gz" -o -iname "*.fna.gz" -o -iname "*.fasta.gz" \\
                -o -iname "*.fq.gz" -o -iname "*.fastq.gz" \\) \\
                -print | sort
        fi
    }

    mapfile -t reference_files < <(resolve_inputs "${reference_dir}")
    mapfile -t read_files < <(resolve_inputs "${reads_dir}")

    if [ "\${#reference_files[@]}" -eq 0 ] && [ "\${#read_files[@]}" -eq 0 ]; then
        echo "No suitable reference or read sequence files were found." >&2
        exit 1
    fi

    plast_args=(
        Build
        -i "graph/${graph_prefix}"
        -k ${params.kmer_length}
        -g ${params.min_length}
        -w ${params.seed_length}
    )

    if [ "\${#reference_files[@]}" -gt 0 ]; then
        plast_args+=( -R "\${reference_files[@]}" )
    fi

    if [ "\${#read_files[@]}" -gt 0 ]; then
        plast_args+=( -S "\${read_files[@]}" )
    fi

    if [ "${params.advanced_index}" = "true" ]; then
        plast_args+=( -a )
    fi

    PLAST "\${plast_args[@]}" \\
        > graph/build.stdout.txt \\
        2> graph/build.stderr.txt

    printf 'PLAST' > graph/build.command.txt

    for argument in "\${plast_args[@]}"; do
        printf ' %q' "\${argument}" >> graph/build.command.txt
    done

    printf '\\n' >> graph/build.command.txt
    printf '%s\\n' "${graph_prefix}" > graph/graph_prefix.txt

    if [ ! -f "graph/${graph_prefix}.gfa" ] && \
       [ ! -f "graph/${graph_prefix}.gfa.gz" ]; then

        echo "Error: Neither '${graph_prefix}.gfa' nor '${graph_prefix}.gfa.gz' exists." >&2
        exit 1
    fi

    if [ ! -f "graph/${graph_prefix}.color.bfg" ]; then
        echo "PLAST Build did not create '${graph_prefix}.color.bfg'." >&2
        exit 1
    fi

    if [ ! -f "graph/${graph_prefix}.idx" ]; then
        echo "PLAST Build did not create '${graph_prefix}.idx'." >&2
        exit 1
    fi
    """
}
