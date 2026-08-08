// Phylogenetic tree construction.
//
// Aligns unique ASV sequences and builds a neighbor-joining tree.
// R uses DECIPHER; Python uses MAFFT + BioPython/scipy.
// Optional — only needed for UniFrac and phylogeny-aware ordinations.

def phyloExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process BUILD_PHYLOGENY {
    tag "phylogeny"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    // Also publish the Newick tree into viz/ as tree.nwk. The viz app fetches
    // data/tree.nwk for its Phylogeny view, and viz/ is what gets deployed as
    // the site's data/ — without this the tree stayed in phylogeny/ and the
    // Phylogeny tab 404'd even when --run_phylogeny was set.
    publishDir "${params.outdir}/viz", mode: 'copy', pattern: 'phylo_tree.nwk',
               saveAs: { 'tree.nwk' }
    storeDir params.store_dir ? "${params.store_dir}/phylogeny" : null

    input:
    path(seqtab)

    output:
    path("phylo_tree.${params.lang == 'python' ? 'nwk' : 'rds'}"), emit: tree
    path("phylo_distances.${phyloExt()}"), emit: distances
    path("phylo_seq_map.${phyloExt()}"), emit: seq_map
    path("phylo_alignment.fasta"), emit: alignment

    script:
    if (params.lang == 'python')
    """
    build_phylogeny.py "${seqtab}" ${task.cpus}
    """
    else
    """
    build_phylogeny.R "${seqtab}" ${task.cpus}
    """
}
