// Gene annotation: Prokka on co-assembly

process PROKKA_ANNOTATE {
    tag "prokka"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-prokka"
    publishDir "${params.outdir}/annotation/prokka", mode: 'copy'

    input:
    path(assembly)

    output:
    path("prokka_out/*.faa"),  emit: proteins
    path("prokka_out/*.gff"),  emit: gff
    path("prokka_out/*.tsv"),  emit: tsv

    script:
    """
    # Stub out tbl2asn -- it hangs on large metagenomes and we don't use its output
    mkdir -p .local/bin
    printf '%s\\n' '#!/bin/sh' 'echo "tbl2asn 25.8"' \\
        'for a in "\$@"; do case "\$p" in -i) f="\$a";; -Z) touch "\$a";; esac; p="\$a"; done' \\
        '[ -n "\$f" ] && touch "\${f%.fsa}.gbf"' 'exit 0' > .local/bin/tbl2asn
    chmod +x .local/bin/tbl2asn
    export PATH="\$PWD/.local/bin:\$PATH"

    prokka \\
        --metagenome \\
        --cpus ${task.cpus} \\
        --evalue 1e-20 \\
        --outdir prokka_out \\
        --force \\
        --quiet \\
        "${assembly}"

    # Clean up unnecessary outputs
    rm -f prokka_out/*.{err,fna,fsa,gbk,log,sqn,txt} 2>/dev/null || true
    rm -f prokka_out/*.tmp.*.faa 2>/dev/null || true

    if [ ! -s prokka_out/PROKKA_*.faa ]; then
        echo "[WARNING] Prokka produced no protein output" >&2
    fi
    """
}
