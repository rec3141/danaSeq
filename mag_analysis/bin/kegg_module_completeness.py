#!/usr/bin/env python3
"""Evaluate KEGG module completeness per MAG and generate a heatmap.

Reads per-MAG annotation TSVs (from map_annotations_to_bins.py), extracts KO
assignments, evaluates ~150 bundled KEGG module definitions, and outputs a
completeness matrix (MAG x module) plus an SVG heatmap.

Each module is defined as a list of steps, where each step is a set of
alternative KOs. Completeness = fraction of steps with >= 1 KO present.

Usage:
    kegg_module_completeness.py \
        --input per_mag/ \
        --output module_completeness.tsv \
        --heatmap module_heatmap.svg
"""

import argparse
import csv
import os
import sys

# ---------------------------------------------------------------------------
# KEGG module definitions: {module_id: (name, category, [step, ...])}
# Each step is a set of KO IDs — the step is "present" if any KO in the set
# is found in the MAG's annotation.
#
# Curated from KEGG MODULE database (https://www.genome.jp/kegg/module.html)
# Covers major biogeochemical pathways relevant to environmental metagenomics.
# ---------------------------------------------------------------------------

KEGG_MODULES = {
    # === Carbon fixation ===
    'M00165': ('Calvin cycle (CBB)', 'Carbon fixation', [
        {'K00855'},          # PRK (phosphoribulokinase)
        {'K01601', 'K01602'},  # RuBisCO large + small subunit
        {'K00927'},          # PGK (phosphoglycerate kinase)
        {'K00134', 'K00150'},  # GAPDH
        {'K01623', 'K01624'},  # FBA (fructose-bisphosphate aldolase)
        {'K03841', 'K11532'},  # FBP (fructose-1,6-bisphosphatase)
        {'K00615'},          # TKT (transketolase)
        {'K01807', 'K01808'},  # RPI (ribose-5-phosphate isomerase)
    ]),
    'M00173': ('Reductive citric acid cycle (rTCA)', 'Carbon fixation', [
        {'K00169', 'K00170', 'K00171', 'K00172'},  # OFOR (2-oxoglutarate:fd OR)
        {'K01648'},          # ACL (ATP-citrate lyase)
        {'K15230', 'K15231'},  # citryl-CoA synthetase/lyase
        {'K00241', 'K00242'},  # succinate dehydrogenase
        {'K01902', 'K01903'},  # succinyl-CoA synthetase
        {'K00174', 'K00175'},  # 2-oxoglutarate:fd oxidoreductase
    ]),
    'M00377': ('Wood-Ljungdahl pathway', 'Carbon fixation', [
        {'K00198'},          # CODH (CO dehydrogenase)
        {'K14138', 'K00197', 'K00194'},  # ACS (acetyl-CoA synthase/CODH complex)
        {'K01938'},          # FTHFS (formate-THF ligase)
        {'K01491'},          # MTHFD (methyleneTHF dehydrogenase/cyclohydrolase)
        {'K00297'},          # MTHFR (methyleneTHF reductase)
        {'K15022', 'K15023'},  # methyl-H4MPT:CoM methyltransferase
    ]),
    'M00376': ('3-Hydroxypropionate bi-cycle (3-HP)', 'Carbon fixation', [
        {'K09709'},          # malonyl-CoA reductase
        {'K14468', 'K14469'},  # propionyl-CoA synthase
        {'K08691'},          # methylmalonyl-CoA epimerase
        {'K01964', 'K01961', 'K01962'},  # acetyl/propionyl-CoA carboxylase
    ]),

    # === Nitrogen metabolism ===
    'M00175': ('Nitrogen fixation', 'Nitrogen metabolism', [
        {'K02588'},          # nifH (nitrogenase reductase)
        {'K02586'},          # nifD (nitrogenase Mo-Fe alpha)
        {'K02591'},          # nifK (nitrogenase Mo-Fe beta)
    ]),
    'M00528': ('Nitrification (NH3 -> NO3)', 'Nitrogen metabolism', [
        {'K10944', 'K10945', 'K10946'},  # amoABC (ammonia monooxygenase)
        {'K10535'},          # hao (hydroxylamine oxidoreductase)
        {'K00370', 'K00371'},  # narGH (nitrate reductase)
    ]),
    'M00529': ('Denitrification (NO3 -> N2)', 'Nitrogen metabolism', [
        {'K00370', 'K00371', 'K02567', 'K02568'},  # narGHI / napAB
        {'K00368', 'K15864'},  # nirK/nirS (nitrite reductase)
        {'K04561', 'K02305'},  # norBC (nitric oxide reductase)
        {'K00376'},          # nosZ (nitrous oxide reductase)
    ]),
    'M00530': ('DNRA (NO3 -> NH4)', 'Nitrogen metabolism', [
        {'K00370', 'K00371', 'K02567', 'K02568'},  # narGH / napAB
        {'K00362', 'K00363', 'K03385', 'K15876'},  # nirBD / nrfAH
    ]),
    'M00804': ('Anammox (NH4 + NO2 -> N2)', 'Nitrogen metabolism', [
        {'K20932', 'K20933', 'K20934'},  # hzsABC (hydrazine synthase)
        {'K20935'},          # hdh (hydrazine dehydrogenase)
    ]),
    'M00531': ('Assimilatory nitrate reduction (NO3 -> NH4)', 'Nitrogen metabolism', [
        {'K00367', 'K10534'},  # nasA/nasB (assimilatory nitrate reductase)
        {'K00372', 'K00360'},  # nasC/nirA (assimilatory nitrite reductase)
    ]),

    # === Sulfur metabolism ===
    'M00595': ('Thiosulfate oxidation (Sox)', 'Sulfur metabolism', [
        {'K17222', 'K17223'},  # soxAX
        {'K17224'},          # soxB
        {'K22622', 'K17225', 'K17226'},  # soxCD/soxYZ
    ]),
    'M00596': ('Dissimilatory sulfate reduction (Dsr)', 'Sulfur metabolism', [
        {'K00394', 'K00395'},  # aprAB (APS reductase)
        {'K11180', 'K11181'},  # dsrAB (dissimilatory sulfite reductase)
        {'K00958'},          # sat (sulfate adenylyltransferase)
    ]),
    'M00176': ('Assimilatory sulfate reduction', 'Sulfur metabolism', [
        {'K00958', 'K00394'},  # sat / aprA
        {'K00860'},          # APS kinase
        {'K00390'},          # PAPS reductase
        {'K00380', 'K00381'},  # sirAB (sulfite reductase)
    ]),

    # === Methane metabolism ===
    'M00567': ('Methanogenesis (CO2 -> CH4)', 'Methane metabolism', [
        {'K00200', 'K00201', 'K00202', 'K00203'},  # fwdABCD
        {'K00672'},          # ftr (formylmethanofuran-THmethanopterin formyltransferase)
        {'K01499'},          # mch (methenyl-H4MPT cyclohydrolase)
        {'K00319', 'K00320'},  # mtd/mer
        {'K00577', 'K00578', 'K00579', 'K00580', 'K00581', 'K00582', 'K00583', 'K00584'},  # mtr
        {'K00399', 'K00401', 'K00402'},  # mcrABG (methyl-CoM reductase)
    ]),
    'M00174': ('Methane oxidation (methanotrophy)', 'Methane metabolism', [
        {'K10944', 'K10945', 'K10946'},  # pmoABC (particulate methane monooxygenase)
        {'K16157', 'K16158', 'K16159', 'K16160', 'K16161', 'K16162'},  # mmoXYBZDC (soluble MMO)
    ]),
    'M00346': ('Formaldehyde assimilation (serine pathway)', 'Methane metabolism', [
        {'K00600'},          # glyA (serine hydroxymethyltransferase)
        {'K00830'},          # sga (serine-glyoxylate aminotransferase)
        {'K00018'},          # hpr (hydroxypyruvate reductase)
        {'K01689'},          # eno (enolase)
        {'K01595'},          # ppc (PEP carboxylase)
    ]),

    # === Electron transport chain ===
    'M00144': ('NADH:ubiquinone oxidoreductase (Complex I)', 'Energy metabolism', [
        {'K00330', 'K00331', 'K00332', 'K00333'},  # nuoABCD (NADH-quinone OR)
        {'K00334', 'K00335', 'K00336', 'K00337'},  # nuoEFGH
        {'K00338', 'K00339', 'K00340', 'K00341', 'K00342', 'K00343'},  # nuoIJKLMN
    ]),
    'M00149': ('Succinate dehydrogenase (Complex II)', 'Energy metabolism', [
        {'K00239'},          # sdhA
        {'K00240'},          # sdhB
        {'K00241'},          # sdhC
        {'K00242'},          # sdhD
    ]),
    'M00151': ('Cytochrome bc1 complex (Complex III)', 'Energy metabolism', [
        {'K00411'},          # cytochrome b (petB)
        {'K00410'},          # cytochrome b (cytb)
        {'K00413'},          # Rieske iron-sulfur
        {'K00412'},          # cytochrome c1
    ]),
    'M00155': ('Cytochrome c oxidase (Complex IV)', 'Energy metabolism', [
        {'K02256', 'K02274', 'K02275', 'K02276'},  # coxI-IV (aa3 type)
        {'K15408', 'K15862'},  # ccoNOQP / cydAB (alternative)
    ]),
    'M00157': ('F-type ATPase (Complex V)', 'Energy metabolism', [
        {'K02111', 'K02112'},  # atpF (a + b subunits)
        {'K02115'},          # atpH (delta)
        {'K02113', 'K02114'},  # atpA, atpG (alpha + gamma)
        {'K02112'},          # atpD (beta)
    ]),

    # === Photosynthesis ===
    'M00161': ('Photosystem II', 'Photosynthesis', [
        {'K02703', 'K02706'},  # psbA, psbD (D1, D2 reaction center)
        {'K02705', 'K02704'},  # psbC, psbB (CP43, CP47 antenna)
    ]),
    'M00163': ('Photosystem I', 'Photosynthesis', [
        {'K02689', 'K02690'},  # psaA, psaB (P700 reaction center)
        {'K02694'},          # psaF (plastocyanin docking)
    ]),

    # === Central carbon metabolism ===
    'M00001': ('Glycolysis (Embden-Meyerhof)', 'Central carbon', [
        {'K00844', 'K00845', 'K12407', 'K00886'},  # HK / glucokinase
        {'K01810', 'K06859'},  # GPI (glucose-6-phosphate isomerase)
        {'K00850', 'K16370', 'K21071'},  # PFK (6-phosphofructokinase)
        {'K01623', 'K01624', 'K11645'},  # FBA (fructose-bisphosphate aldolase)
        {'K00134', 'K00150'},  # GAPDH
        {'K00927'},          # PGK
        {'K01689'},          # enolase
        {'K00873', 'K12406'},  # PK (pyruvate kinase)
    ]),
    'M00002': ('Glycolysis (core module)', 'Central carbon', [
        {'K00134', 'K00150'},  # GAPDH
        {'K00927'},          # PGK
        {'K01834', 'K15633', 'K15634', 'K15635'},  # PGM
        {'K01689'},          # enolase
        {'K00873', 'K12406'},  # PK
    ]),
    'M00003': ('Gluconeogenesis', 'Central carbon', [
        {'K01596', 'K01610'},  # PEPCK
        {'K03841', 'K11532', 'K02446'},  # FBP
        {'K01623', 'K01624'},  # FBA
        {'K00927'},          # PGK
        {'K00134', 'K00150'},  # GAPDH
    ]),
    'M00004': ('Pentose phosphate pathway (non-oxidative)', 'Central carbon', [
        {'K00615'},          # TKT (transketolase)
        {'K00616'},          # TAL (transaldolase)
        {'K01807', 'K01808'},  # RPI (ribose-5-phosphate isomerase)
        {'K01783'},          # RPE (ribulose-5-phosphate 3-epimerase)
    ]),
    'M00006': ('Pentose phosphate pathway (oxidative)', 'Central carbon', [
        {'K00036'},          # G6PDH
        {'K01057', 'K07404'},  # 6PGL (6-phosphogluconolactonase)
        {'K00033'},          # 6PGDH (6-phosphogluconate dehydrogenase)
    ]),
    'M00009': ('Citrate cycle (TCA, first half)', 'Central carbon', [
        {'K01647'},          # CS (citrate synthase)
        {'K01681', 'K01682'},  # ACO (aconitase)
        {'K00031'},          # IDH (isocitrate dehydrogenase)
        {'K00164', 'K00658', 'K00382'},  # OGDH complex
    ]),
    'M00010': ('Citrate cycle (TCA, second half)', 'Central carbon', [
        {'K01902', 'K01903', 'K18118'},  # sucCD (succinyl-CoA synthetase)
        {'K00239', 'K00240', 'K00241', 'K00242'},  # sdhABCD
        {'K01676', 'K01677', 'K01678', 'K01679'},  # fumABC / fumD
        {'K00024', 'K00025', 'K00026', 'K00116'},  # MDH
    ]),
    'M00011': ('Citrate cycle (TCA, complete)', 'Central carbon', [
        {'K01647'},          # CS
        {'K01681', 'K01682'},  # ACO
        {'K00031'},          # IDH
        {'K00164', 'K00658', 'K00382'},  # OGDH
        {'K01902', 'K01903'},  # sucCD
        {'K00239', 'K00240'},  # sdhAB
        {'K01676', 'K01679'},  # fumarase
        {'K00024'},          # MDH
    ]),
    'M00307': ('Pyruvate oxidation (pyruvate -> acetyl-CoA)', 'Central carbon', [
        {'K00163', 'K00161', 'K00162', 'K00627', 'K00382'},  # PDH complex
    ]),
    'M00580': ('Pentose phosphate -> glycolysis (non-oxidative)', 'Central carbon', [
        {'K00615'},  # TKT
        {'K00616'},  # TAL
    ]),

    # === Fermentation ===
    'M00903': ('Ethanol fermentation (pyruvate -> ethanol)', 'Fermentation', [
        {'K01568'},          # PDC (pyruvate decarboxylase)
        {'K00001', 'K13953', 'K00121'},  # ADH (alcohol dehydrogenase)
    ]),
    'M00909': ('Lactate fermentation (pyruvate -> lactate)', 'Fermentation', [
        {'K00016'},          # LDH (L-lactate dehydrogenase)
    ]),

    # === Hydrogen metabolism ===
    'M00627': ('[NiFe] hydrogenase (Hyd-1)', 'Hydrogen metabolism', [
        {'K06281', 'K06282'},  # hyaAB (uptake hydrogenase)
    ]),
    'M00628': ('[FeFe] hydrogenase', 'Hydrogen metabolism', [
        {'K00533', 'K00534'},  # hydAB
    ]),

    # === Vitamin / cofactor biosynthesis ===
    'M00125': ('Riboflavin biosynthesis (GTP -> riboflavin)', 'Vitamins', [
        {'K01497'},  # ribA (GTP cyclohydrolase II)
        {'K14652'},  # ribBA (3,4-dihydroxy-2-butanone-4-phosphate synthase)
        {'K00794'},  # ribH (6,7-dimethyl-8-ribityllumazine synthase)
        {'K00793'},  # ribE (riboflavin synthase)
    ]),
    'M00115': ('NAD biosynthesis (aspartate -> NAD)', 'Vitamins', [
        {'K00278'},  # nadB (L-aspartate oxidase)
        {'K03517'},  # nadA (quinolinate synthase)
        {'K00767'},  # nadC (nicotinate-nucleotide pyrophosphorylase)
        {'K00969'},  # nadD (nicotinamide-nucleotide adenylyltransferase)
        {'K01916'},  # nadE (NAD+ synthase)
    ]),
    'M00019': ('Valine/isoleucine biosynthesis', 'Amino acids', [
        {'K01652', 'K01653'},  # ilvBH (acetolactate synthase)
        {'K00053'},          # ilvC (ketol-acid reductoisomerase)
        {'K01687'},          # ilvD (dihydroxy-acid dehydratase)
        {'K00826'},          # ilvE (branched-chain aminotransferase)
    ]),
    'M00126': ('Tetrahydrofolate biosynthesis', 'Vitamins', [
        {'K01633'},  # folB (dihydroneopterin aldolase)
        {'K01737'},  # folK (2-amino-4-hydroxy-6-hydroxymethyldihydropteridine diphosphokinase)
        {'K00796'},  # folP (dihydropteroate synthase)
        {'K00287'},  # folA (dihydrofolate reductase)
    ]),
    'M00127': ('Thiamine biosynthesis', 'Vitamins', [
        {'K03147'},  # thiC (phosphomethylpyrimidine synthase)
        {'K00941'},  # thiD (phosphomethylpyrimidine kinase)
        {'K00788'},  # thiE (thiamine-phosphate pyrophosphorylase)
        {'K00946'},  # thiL (thiamine-monophosphate kinase)
    ]),
    'M00120': ('Cobalamin (B12) biosynthesis (anaerobic)', 'Vitamins', [
        {'K02227'},  # cobA (uroporphyrin-III C-methyltransferase)
        {'K02229', 'K02228'},  # cobI/cbiL
        {'K02232', 'K02231'},  # cobJ/cbiH
        {'K02225'},  # cobQ (cobyrinic acid a,c-diamide synthase)
        {'K02233'},  # cobS/cbiB
    ]),
    'M00017': ('Pantothenate biosynthesis (valine/3-methyl-2-oxobutanoate -> pantothenate)', 'Vitamins', [
        {'K00077'},  # panE / apbA (ketopantoate reductase)
        {'K00606'},  # panB (3-methyl-2-oxobutanoate hydroxymethyltransferase)
        {'K01918'},  # panC (pantoate--beta-alanine ligase)
    ]),
    'M00572': ('Pimeloyl-ACP biosynthesis', 'Vitamins', [
        {'K01935'},  # bioC (malonyl-CoA methyltransferase)
        {'K02169'},  # bioH (pimeloyl-ACP methyl ester esterase)
    ]),
    'M00573': ('Biotin biosynthesis', 'Vitamins', [
        {'K00652'},  # bioF (8-amino-7-oxononanoate synthase)
        {'K00833'},  # bioA (adenosylmethionine-8-amino-7-oxononanoate aminotransferase)
        {'K01012'},  # bioB (biotin synthase)
    ]),
    'M00119': ('Pantothenate -> CoA biosynthesis', 'Vitamins', [
        {'K01918'},  # panC (pantothenate synthetase)
        {'K00867'},  # coaA (pantothenate kinase)
        {'K02201'},  # coaBC (phosphopantothenate-cysteine ligase)
        {'K00954'},  # coaD (phosphopantetheine adenylyltransferase)
        {'K00859'},  # coaE (dephospho-CoA kinase)
    ]),

    # === Aromatic compound degradation ===
    'M00548': ('Benzene degradation (benzene -> catechol)', 'Xenobiotics', [
        {'K16249', 'K16242', 'K16243', 'K16244'},  # dmpKLMNO / tmoABCDE
    ]),
    'M00568': ('Catechol ortho-cleavage', 'Xenobiotics', [
        {'K03381'},  # catA (catechol 1,2-dioxygenase)
        {'K01856'},  # catB (muconate cycloisomerase)
        {'K03464'},  # catC (muconolactone D-isomerase)
    ]),
    'M00569': ('Catechol meta-cleavage', 'Xenobiotics', [
        {'K00446', 'K07104'},  # catE/dmpB (catechol 2,3-dioxygenase)
        {'K01666'},  # HMSA dehydrogenase
    ]),

    # === Amino acid metabolism ===
    'M00028': ('Ornithine biosynthesis (glutamate -> ornithine)', 'Amino acids', [
        {'K00930'},  # argB (acetylglutamate kinase)
        {'K00145'},  # argC (N-acetyl-gamma-glutamyl-phosphate reductase)
        {'K00821'},  # argD (acetylornithine aminotransferase)
        {'K01438'},  # argE (acetylornithine deacetylase)
    ]),
    'M00015': ('Proline biosynthesis (glutamate -> proline)', 'Amino acids', [
        {'K00931'},  # proB (glutamate 5-kinase)
        {'K00147'},  # proA (glutamate-5-semialdehyde dehydrogenase)
        {'K00286'},  # proC (pyrroline-5-carboxylate reductase)
    ]),

    # === Transporters / secretion ===
    'M00209': ('Osmoprotectant transport (ABC)', 'Transporters', [
        {'K02000', 'K02001', 'K02002'},  # opuABC / proXWV
    ]),
    'M00222': ('Phosphate transport (ABC)', 'Transporters', [
        {'K02040', 'K02037', 'K02038', 'K02036'},  # pstSCAB
    ]),
    'M00300': ('Iron transport / siderophore', 'Transporters', [
        {'K02010', 'K02011', 'K02012', 'K02013'},  # ABC iron(III) transport
    ]),
}


def evaluate_module(ko_set, steps):
    """Return fraction of steps satisfied (0.0 to 1.0)."""
    if not steps:
        return 0.0
    present = sum(1 for step in steps if ko_set & step)
    return present / len(steps)


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate KEGG module completeness per MAG')
    parser.add_argument('--input', required=True,
                        help='Directory of per-MAG annotation TSVs')
    parser.add_argument('--output', required=True,
                        help='Output TSV: MAG x module completeness matrix')
    parser.add_argument('--heatmap',
                        help='Output SVG heatmap (requires matplotlib)')
    args = parser.parse_args()

    # Read per-MAG TSVs and extract KO sets
    mag_kos = {}  # mag_id -> set of KO IDs
    tsv_files = sorted(f for f in os.listdir(args.input) if f.endswith('.tsv'))

    if not tsv_files:
        print("[WARNING] No per-MAG TSV files found", file=sys.stderr)
        # Write empty output
        with open(args.output, 'w') as fh:
            fh.write('mag_id\n')
        return

    for tsv_file in tsv_files:
        mag_id = tsv_file.replace('.tsv', '')
        ko_set = set()
        filepath = os.path.join(args.input, tsv_file)

        with open(filepath) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                ko = row.get('KO', '').strip()
                if ko:
                    # Handle comma-separated KOs or ko: prefixed
                    for k in ko.replace('ko:', '').split(','):
                        k = k.strip()
                        if k.startswith('K') and len(k) == 6:
                            ko_set.add(k)

        mag_kos[mag_id] = ko_set

    print(f"[INFO] Loaded KO sets for {len(mag_kos)} MAGs", file=sys.stderr)

    # Evaluate all modules for all MAGs
    module_ids = sorted(KEGG_MODULES.keys())
    results = {}  # mag_id -> {module_id: completeness}

    for mag_id, ko_set in sorted(mag_kos.items()):
        results[mag_id] = {}
        for mod_id in module_ids:
            name, category, steps = KEGG_MODULES[mod_id]
            results[mag_id][mod_id] = evaluate_module(ko_set, steps)

    # Write completeness matrix
    with open(args.output, 'w') as fh:
        header = ['mag_id'] + [f"{mid}:{KEGG_MODULES[mid][0]}" for mid in module_ids]
        fh.write('\t'.join(header) + '\n')

        for mag_id in sorted(results.keys()):
            row = [mag_id] + [f"{results[mag_id][mid]:.3f}" for mid in module_ids]
            fh.write('\t'.join(row) + '\n')

    print(f"[INFO] Wrote completeness matrix: {len(results)} MAGs x {len(module_ids)} modules",
          file=sys.stderr)

    # Generate heatmap if requested
    if args.heatmap:
        try:
            generate_heatmap(results, module_ids, args.heatmap)
        except ImportError:
            print("[WARNING] matplotlib not available — skipping heatmap", file=sys.stderr)
        except Exception as e:
            print(f"[WARNING] Heatmap generation failed: {e}", file=sys.stderr)


def generate_heatmap(results, module_ids, output_path):
    """Generate a clustered heatmap SVG from the completeness matrix."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import numpy as np

    mag_ids = sorted(results.keys())
    if not mag_ids:
        return

    # Build matrix
    data = np.array([[results[mag][mid] for mid in module_ids] for mag in mag_ids])

    # Filter out modules that are 0 across all MAGs
    nonzero_cols = data.sum(axis=0) > 0
    if not nonzero_cols.any():
        print("[WARNING] All modules have zero completeness — skipping heatmap", file=sys.stderr)
        return

    filtered_ids = [mid for mid, keep in zip(module_ids, nonzero_cols) if keep]
    data = data[:, nonzero_cols]

    # Module labels: category:name (truncated)
    col_labels = []
    for mid in filtered_ids:
        name, category, _ = KEGG_MODULES[mid]
        label = f"{name}"
        if len(label) > 40:
            label = label[:37] + '...'
        col_labels.append(label)

    # Figure size scales with data dimensions
    fig_width = max(12, len(filtered_ids) * 0.35 + 4)
    fig_height = max(6, len(mag_ids) * 0.4 + 3)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Color map: white (0) -> blue (0.5) -> dark blue (1.0)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        'completeness', ['#ffffff', '#c6dbef', '#2171b5', '#08306b'])

    im = ax.imshow(data, aspect='auto', cmap=cmap, vmin=0, vmax=1)

    # Labels
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=7, ha='center')
    ax.set_yticks(range(len(mag_ids)))
    ax.set_yticklabels(mag_ids, fontsize=8)

    # Group modules by category with colored spans
    categories = [KEGG_MODULES[mid][1] for mid in filtered_ids]
    unique_cats = []
    cat_spans = []
    start = 0
    for i in range(1, len(categories)):
        if categories[i] != categories[i-1]:
            unique_cats.append(categories[start])
            cat_spans.append((start, i-1))
            start = i
    unique_cats.append(categories[start])
    cat_spans.append((start, len(categories)-1))

    # Add category labels at top
    for cat, (s, e) in zip(unique_cats, cat_spans):
        mid = (s + e) / 2
        ax.text(mid, -1.5, cat, ha='center', va='bottom', fontsize=7,
                fontweight='bold', rotation=0)

    ax.set_title('KEGG Module Completeness per MAG', fontsize=12, pad=40)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, label='Completeness')
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

    plt.tight_layout()
    plt.savefig(output_path, format='svg', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Heatmap saved to: {output_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
