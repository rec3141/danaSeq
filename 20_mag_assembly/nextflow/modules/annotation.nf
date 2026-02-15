// Gene annotation: Prokka on co-assembly

process PROKKA_ANNOTATE {
    tag "prokka"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-prokka"
    publishDir "${params.outdir}/annotation/prokka", mode: 'copy',
        saveAs: { fn -> fn.minus('prokka_out/') }

    input:
    path(assembly)

    output:
    path("prokka_out/*.faa"),  emit: proteins
    path("prokka_out/*.gff"),  emit: gff
    path("prokka_out/*.gbk"),  emit: gbk
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
    rm -f prokka_out/*.{err,fna,fsa,log,sqn,txt} 2>/dev/null || true
    rm -f prokka_out/*.tmp.*.faa 2>/dev/null || true

    if [ ! -s prokka_out/PROKKA_*.faa ]; then
        echo "[WARNING] Prokka produced no protein output" >&2
    fi
    """
}

// Gene annotation: Bakta CDS-only (fast path — minutes on large metagenomes)
// Skips ncRNA/tRNA/CRISPR/sORF scanning (cmscan is the bottleneck on large assemblies).
// Produces .faa + .gff3 for all downstream tools (Kaiju, DefenseFinder, metabolism, etc.)

process BAKTA_CDS {
    tag "bakta_cds"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-bakta"
    publishDir "${params.outdir}/annotation/bakta/cds", mode: 'copy',
        saveAs: { fn -> fn.minus('bakta_out/') }

    input:
    path(assembly)

    output:
    path("bakta_out/annotation.faa"),  emit: proteins
    path("bakta_out/annotation.gff3"), emit: gff
    path("bakta_out/annotation.tsv"),  emit: tsv
    path("bakta_out/annotation.hypotheticals.faa"), emit: hypotheticals, optional: true

    script:
    """
    bakta \
        --db "${params.bakta_db}" \
        --meta \
        --threads ${task.cpus} \
        --output bakta_out \
        --prefix annotation \
        --force \
        --skip-trna --skip-tmrna --skip-rrna \
        --skip-ncrna --skip-ncrna-region \
        --skip-crispr --skip-sorf \
        --skip-gap --skip-ori --skip-plot \
        "${assembly}"

    if [ ! -s bakta_out/annotation.faa ]; then
        echo "[WARNING] Bakta CDS produced no protein output" >&2
    fi
    """
}

// Gene annotation: Bakta full (slow path — hours on large metagenomes)
// Complete annotation including ncRNA, tRNA, CRISPR, sORF, etc.
// Runs in parallel with downstream tools; does not block the pipeline.

process BAKTA_FULL {
    tag "bakta_full"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-bakta"
    publishDir "${params.outdir}/annotation/bakta/full", mode: 'copy',
        saveAs: { fn -> fn.minus('bakta_out/') }

    input:
    path(assembly)

    output:
    path("bakta_out/annotation.faa"),  emit: proteins
    path("bakta_out/annotation.gff3"), emit: gff
    path("bakta_out/annotation.tsv"),  emit: tsv
    path("bakta_out/annotation.hypotheticals.faa"), emit: hypotheticals, optional: true
    path("bakta_out/annotation.gbff"), emit: gbff, optional: true
    path("bakta_out/annotation.png"),  emit: plot, optional: true

    script:
    """
    bakta \
        --db "${params.bakta_db}" \
        --meta \
        --threads ${task.cpus} \
        --output bakta_out \
        --prefix annotation \
        --force \
        "${assembly}"

    if [ ! -s bakta_out/annotation.faa ]; then
        echo "[WARNING] Bakta full produced no protein output" >&2
    fi
    """
}
