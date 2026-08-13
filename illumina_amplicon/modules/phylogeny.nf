// Phylogenetic tree construction.
//
// Aligns unique ASV sequences and builds a neighbor-joining tree.
// R uses DECIPHER; Python uses MAFFT + BioPython/scipy.
// Optional — only needed for UniFrac and phylogeny-aware ordinations.


process BUILD_PHYLOGENY {
    tag "phylogeny"
    label 'process_high'
    conda "${projectDir}/envs/python.yml"
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
    path("phylo_tree.${'nwk'}"), emit: tree
    path("phylo_distances.pkl"), emit: distances
    path("phylo_seq_map.pkl"), emit: seq_map
    path("phylo_alignment.fasta"), emit: alignment

    script:
    """
    build_phylogeny.py "${seqtab}" ${task.cpus}
    """
}
