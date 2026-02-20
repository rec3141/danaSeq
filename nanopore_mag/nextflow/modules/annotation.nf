// Gene annotation: Prokka on co-assembly

process PROKKA_ANNOTATE {
    tag "prokka"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-prokka"
    publishDir "${params.outdir}/annotation/prokka", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/annotation/prokka" : null

    input:
    path(assembly)

    output:
    path("annotation.faa"),  emit: proteins
    path("annotation.gff"),  emit: gff
    path("annotation.gbk"),  emit: gbk
    path("annotation.tsv"),  emit: tsv

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
        --prefix annotation \\
        --outdir prokka_out \\
        --force \\
        --quiet \\
        "${assembly}"

    # Clean up unnecessary outputs
    rm -f prokka_out/*.{err,fna,fsa,log,sqn,txt} 2>/dev/null || true
    rm -f prokka_out/*.tmp.*.faa 2>/dev/null || true

    if [ ! -s prokka_out/annotation.faa ]; then
        echo "[WARNING] Prokka produced no protein output" >&2
    fi

    cp prokka_out/annotation.faa prokka_out/annotation.gff prokka_out/annotation.gbk prokka_out/annotation.tsv .
    """
}

// Gene annotation: Bakta basic (fast path — minutes on large metagenomes)
// Uses the Bakta light database (auto-derived from bakta_db path).
// Skips ncRNA/tRNA/CRISPR/sORF scanning (cmscan is the bottleneck on large assemblies).
// Produces .faa + .gff3 for all downstream tools (Kaiju, DefenseFinder, metabolism, etc.)

process BAKTA_BASIC {
    tag "bakta_basic"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-bakta"
    publishDir "${params.outdir}/annotation/bakta/basic", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/annotation/bakta/basic" : null

    input:
    path(assembly)

    output:
    path("annotation.faa"),  emit: proteins
    path("annotation.gff3"), emit: gff
    path("annotation.tsv"),  emit: tsv
    path("annotation.hypotheticals.faa"), emit: hypotheticals, optional: true

    script:
    // Derive light DB: bakta_db points to full DB (e.g. .../bakta/full/db),
    // light DB is its sibling (e.g. .../bakta/light/db-light)
    def light_db = file(params.bakta_db).parent.parent.resolve('light/db-light')
    """
    if [ ! -d "${light_db}" ]; then
        echo "[ERROR] Bakta light database not found at ${light_db}" >&2
        echo "[INFO]  Expected layout: <bakta_root>/full/db and <bakta_root>/light/db-light" >&2
        echo "[INFO]  Download with: bakta_db download --output \$(dirname \$(dirname ${params.bakta_db})) --type light" >&2
        exit 1
    fi

    bakta \
        --db "${light_db}" \
        --meta \
        --keep-contig-headers \
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
        echo "[WARNING] Bakta basic produced no protein output" >&2
    fi

    cp bakta_out/annotation.faa bakta_out/annotation.tsv .
    # Strip ##FASTA section from GFF3 (embedded sequences bloat the file and confuse parsers)
    sed '/^##FASTA/Q' bakta_out/annotation.gff3 > annotation.gff3
    cp bakta_out/annotation.hypotheticals.faa . 2>/dev/null || true
    """
}

// Gene annotation: Bakta extra (full DB — CDS + CRISPR + sORF + expert systems)
// Uses the Bakta full database. Skips RNA annotation (handled by RNA_CLASSIFY module).
// Runs in parallel with downstream tools; does not block the pipeline.

process BAKTA_EXTRA {
    tag "bakta_extra"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-bakta"
    publishDir "${params.outdir}/annotation/bakta/extra", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/annotation/bakta/extra" : null

    input:
    path(assembly)

    output:
    path("annotation.faa"),  emit: proteins
    path("annotation.gff3"), emit: gff
    path("annotation.tsv"),  emit: tsv
    path("annotation.hypotheticals.faa"), emit: hypotheticals, optional: true
    path("annotation.gbff"), emit: gbff, optional: true
    path("annotation.png"),  emit: plot, optional: true

    script:
    """
    bakta \
        --db "${params.bakta_db}" \
        --meta \
        --keep-contig-headers \
        --threads ${task.cpus} \
        --output bakta_out \
        --prefix annotation \
        --force \
        --skip-trna --skip-tmrna --skip-rrna \
        --skip-ncrna --skip-ncrna-region \
        "${assembly}"

    if [ ! -s bakta_out/annotation.faa ]; then
        echo "[WARNING] Bakta full produced no protein output" >&2
    fi

    cp bakta_out/annotation.faa bakta_out/annotation.tsv .
    # Strip ##FASTA section from GFF3 (embedded sequences bloat the file and confuse parsers)
    sed '/^##FASTA/Q' bakta_out/annotation.gff3 > annotation.gff3
    cp bakta_out/annotation.hypotheticals.faa bakta_out/annotation.gbff bakta_out/annotation.png . 2>/dev/null || true
    """
}
