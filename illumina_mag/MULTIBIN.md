 Mapping Strategies

  1. Cross-mapping to per-sample assemblies (what your original METTA 51_/52_/56_ scripts do) — Map reads from all N samples to each sample's assembly. Produces an N-column depth matrix per assembly. Best bin recovery for individual
  assemblies.
  2. Co-assembly + cross-mapping — Pool all reads into one assembly, map all samples back. Simpler DAG but loses rare/sample-specific organisms. The plan already has --coassembly placeholder.
  3. Hybrid — Per-sample assemblies, then cross-map a subset (e.g., samples from the same site/cruise) rather than all-vs-all.

  Binning Tools

  ┌──────────┬───────────────────────────────┬────────────────────────────────┬──────────────────────────────────────────────────────────┐
  │   Tool   │           Approach            │         Multi-sample?          │                          Notes                           │
  ├──────────┼───────────────────────────────┼────────────────────────────────┼──────────────────────────────────────────────────────────┤
  │ MetaBAT2 │ Composition + coverage        │ Yes (multi-BAM depth matrix)   │ Already in pipeline. Workhorse.                          │
  ├──────────┼───────────────────────────────┼────────────────────────────────┼──────────────────────────────────────────────────────────┤
  │ CONCOCT  │ Composition + coverage        │ Yes (designed for it)          │ Good complement to MetaBAT2                              │
  ├──────────┼───────────────────────────────┼────────────────────────────────┼──────────────────────────────────────────────────────────┤
  │ MaxBin2  │ Composition + coverage        │ Yes (multiple abundance files) │ Older but still used                                     │
  ├──────────┼───────────────────────────────┼────────────────────────────────┼──────────────────────────────────────────────────────────┤
  │ VAMB     │ Variational autoencoder       │ Yes (purpose-built)            │ State-of-the-art for multi-sample, needs GPU or patience │
  ├──────────┼───────────────────────────────┼────────────────────────────────┼──────────────────────────────────────────────────────────┤
  │ SemiBin2 │ Self-supervised deep learning │ Yes (multi-sample mode)        │ Very strong recent performer                             │
  └──────────┴───────────────────────────────┴────────────────────────────────┴──────────────────────────────────────────────────────────┘

  Bin Refinement / Consensus

  ┌──────────┬──────────────────────────────────────────────────────────────────────────────────┐
  │   Tool   │                                     Purpose                                      │
  ├──────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ DAS Tool │ Picks best bins from multiple binners (MetaBAT2 + CONCOCT + MaxBin2 → consensus) │
  ├──────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ CheckM2  │ Quality assessment (completeness/contamination) for filtering                    │
  ├──────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ dRep     │ Dereplicate bins across samples at species-level (ANI threshold)                 │
  ├──────────┼──────────────────────────────────────────────────────────────────────────────────┤
  │ GTDB-Tk  │ Taxonomic classification of final MAGs                                           │
  └──────────┴──────────────────────────────────────────────────────────────────────────────────┘

  Common Production Workflows

  Minimal (what to add next): Cross-mapping + MetaBAT2 with multi-sample depths — matches your original bash pipeline.

  Standard: Cross-mapping + MetaBAT2 + CONCOCT + MaxBin2 → DAS Tool consensus → CheckM2 filtering

  State-of-the-art: Cross-mapping + VAMB or SemiBin2 (often outperform the 3-binner+DAS Tool approach alone) → CheckM2 → dRep → GTDB-Tk

  The most impactful next step would be implementing the cross-mapping module, since that unlocks multi-sample depth for any binner and matches what your original METTA scripts already do. Which direction interests you?
