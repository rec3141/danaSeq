"""MetaCyc pathway predictions for mag_analysis (entry point: bin/metacyc_scorer.py).

ptdat       Pathway Tools attribute-value parser and the MetaCyc pathway/reaction model
evidence    protein -> Rhea -> MetaCyc reaction assignments, chunked by protein
pathways    PathoLogic's published pathway rules and MinPath-style parsimony
confidence  evidence-weighted confidence and tiers for each call
taxonomy    GTDB lineage -> NCBI taxids for MetaCyc taxonomic ranges
refs        builder for the open references (Rhea, UniProt DIAMOND reference, GTDB)
install     installs the user's licensed MetaCyc from flat files, an installer or a .sif
"""
