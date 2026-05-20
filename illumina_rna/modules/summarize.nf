// Cross-sample summaries: merged gene-count matrix, mapping-rate matrix,
// and the viz JSONs.

process MERGE_GENE_COUNTS {
    tag "${ref.name}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-rna-rnaseq"
    publishDir "${params.outdir}/expression/${ref.name}", mode: 'copy'

    input:
    tuple val(ref), path(per_sample_counts, stageAs: 'counts_*.tsv')

    output:
    tuple val(ref), path("${ref.name}.gene_counts.tsv"), emit: counts

    script:
    """
    python3 - <<PY
import glob, os, pandas as pd
files = sorted(glob.glob("counts_*.tsv"))
dfs = []
for f in files:
    df = pd.read_csv(f, sep="\\t", comment="#")
    if df.empty or len(df.columns) < 7:
        continue
    sample_col = df.columns[-1]
    sub = df[["Geneid", sample_col]].set_index("Geneid")
    dfs.append(sub)
if not dfs:
    pd.DataFrame(columns=["Geneid"]).to_csv("${ref.name}.gene_counts.tsv", sep="\\t", index=False)
else:
    merged = pd.concat(dfs, axis=1).fillna(0).astype(int)
    merged.to_csv("${ref.name}.gene_counts.tsv", sep="\\t")
PY
    """
}

process VIZ_PREPROCESS {
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-rna-rnaseq"
    publishDir "${params.outdir}/viz", mode: 'copy'

    input:
    path(idxstats_files,  stageAs: 'idxstats/*')
    path(flagstat_files,  stageAs: 'flagstat/*')
    path(covstats_files,  stageAs: 'covstats/*')
    path(read_counts,     stageAs: 'readcounts/*')
    path(gene_counts,     stageAs: 'genecounts/*')

    output:
    path("overview.json"),       emit: overview
    path("samples.json"),        emit: samples
    path("references.json"),     emit: references
    path("expression.json.gz"),  emit: expression
    path("read_flow.json"),      emit: read_flow

    script:
    """
    python3 "${projectDir}/bin/build_viz_jsons.py" \\
        --idxstats-dir idxstats \\
        --flagstat-dir flagstat \\
        --covstats-dir covstats \\
        --readcounts-dir readcounts \\
        --genecounts-dir genecounts \\
        --out .
    """
}
