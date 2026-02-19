// Prokka gene annotation in metagenome mode

process PROKKA_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-prokka"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/prokka", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}/PROKKA_*.faa"),  emit: proteins
    tuple val(meta), path("${meta.id}/*.tsv"),  emit: tsv
    tuple val(meta), path("${meta.id}/*.gff"),  emit: gff
    tuple val(meta), path("${meta.id}"),        emit: prokka_dir

    script:
    """
    # Stub out tbl2asn -- it hangs on large metagenomes and we don't use its output
    # Prokka checks version via: tbl2asn - | grep '^tbl2asn' then matches /tbl2asn\\s+(\\d+\\.\\d+)/
    # Prokka also runs tbl2asn for annotation and expects .gbf output file (sed .gbf -> .gbk)
    mkdir -p .local/bin
    printf '%s\\n' '#!/bin/sh' 'echo "tbl2asn 25.8"' 'for a in "\$@"; do case "\$p" in -i) f="\$a";; -Z) touch "\$a";; esac; p="\$a"; done' '[ -n "\$f" ] && touch "\${f%.fsa}.gbf"' 'exit 0' > .local/bin/tbl2asn
    chmod +x .local/bin/tbl2asn
    export PATH="\$PWD/.local/bin:\$PATH"

    prokka \
        --metagenome \
        --fast \
        --cpus 1 \
        --evalue 1e-20 \
        --outdir "${meta.id}" \
        --force \
        --quiet \
        "${fasta}"

    # Clean up unnecessary Prokka outputs to save space
    rm -f "${meta.id}"/*.{err,fna,fsa,gbk,log,sqn,txt} 2>/dev/null || true
    # Remove temp files that can interfere with downstream glob matching
    rm -f "${meta.id}"/*.tmp.*.faa 2>/dev/null || true
    """
}
