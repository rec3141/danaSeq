#!/usr/bin/env bash
################################################################################
#                                                                              #
#  🦠  AGENT: THE MICROBIAL ECOLOGIST  🦠                                      #
#                                                                              #
#  "Everything you think you know about ecology? It started with microbes."   #
#                                                                              #
################################################################################

readonly MAGENTA='\033[0;35m'
readonly GREEN='\033[0;32m'
readonly BOLD='\033[1m'
readonly NC='\033[0m'

cat << "EOF"

     ╔═══════════════════════════════════════════════════╗
     ║  🦠    THE MICROBIAL ECOLOGIST SPEAKS    🦠      ║
     ╚═══════════════════════════════════════════════════╝

     "I study the unseen majority. They're small,
      but they run the planet. Let me show you how."

     ╔═══════════════════════════════════════════════════╗

EOF

echo -e "${MAGENTA}${BOLD}WHO AM I?${NC}\n"
echo "Dr. Petra Microbe, Microbial Ecologist"
echo "Specialties: Community ecology, metabolic modeling, evolution"
echo "Cultures isolated: 247 (but 10,000+ remain unculturable)"
echo ""

echo -e "${MAGENTA}${BOLD}WHAT I CARE ABOUT:${NC}\n"

echo -e "${GREEN}🔬 Community Structure${NC}"
echo "   • Who's there? (taxonomy, diversity)"
echo "   • How many? (abundance, dominance)"
echo "   • Who's with whom? (co-occurrence, networks)"
echo "   • Who's active? (transcriptomics, not just DNA!)"
echo ""

echo -e "${GREEN}⚡ Metabolic Interactions${NC}"
echo "   • Syntrophy (metabolic handshakes)"
echo "   • Competition (resource limitation)"
echo "   • Predation (viruses, protists)"
echo "   • Cooperation (public goods production)"
echo ""

echo -e "${GREEN}🌍 Ecosystem Function${NC}"
echo "   • What are they doing? (pathways, reactions)"
echo "   • At what rate? (fluxes, turnover)"
echo "   • Why here, not there? (niche partitioning)"
echo "   • How do they respond? (perturbations, succession)"
echo ""

echo -e "${MAGENTA}${BOLD}KEY ECOLOGICAL CONCEPTS FOR YOUR DATA:${NC}\n"

echo -e "${GREEN}1. THE RARE BIOSPHERE${NC}"
echo "   • Most taxa are <0.1% abundance"
echo "   • High diversity, low abundance"
echo "   • Seed bank for future blooms"
echo "   • Genetic reservoir"
echo ""
echo "   Your MAGs will miss them (need deeper sequencing)!"
echo "   But they matter: rare today = dominant tomorrow"
echo ""

echo -e "${GREEN}2. NICHE PARTITIONING${NC}"
echo "   • Every microbe has a niche"
echo "   • Similar niches = competition = exclusion"
echo "   • Hutchinsonian paradox: Why so many species?"
echo ""
echo "   Example: Prochlorococcus ecotypes"
echo "   ├─ High-light ecotype (surface)"
echo "   ├─ Low-light ecotype (DCM)"
echo "   └─ Same genus, different niches, coexist!"
echo ""

echo -e "${GREEN}3. BLOOM DYNAMICS${NC}"
echo "   • Disturbance → Resource pulse → Boom → Bust"
echo "   • r-selected: Fast growth, low efficiency"
echo "   • K-selected: Slow growth, high efficiency"
echo ""
echo "   Watch for:"
echo "   ├─ Cyanobacteria blooms (eutrophication)"
echo "   ├─ Flavobacteria blooms (phytoplankton die-off)"
echo "   └─ Vibrio blooms (warm + nutrients)"
echo ""

echo -e "${GREEN}4. THE BLACK QUEEN HYPOTHESIS${NC}"
echo "   • Genes cost energy"
echo "   • Lose non-essential genes = grow faster"
echo "   • Depend on neighbors for 'public goods'"
echo ""
echo "   Example: SAR11"
echo "   ├─ Lost biotin synthesis genes"
echo "   ├─ Depends on others making biotin"
echo "   └─ Streamlined genome = competitive advantage"
echo ""

echo -e "${GREEN}5. MICROBIAL LOOP${NC}"
echo ""
cat << "EOF"
   Phytoplankton → Dissolved Organic Matter (DOM)
                           ↓
                    Bacteria consume DOM
                           ↓
                   Protists eat bacteria
                           ↓
                Zooplankton eat protists
                           ↓
                    Fish eat zooplankton

   Without bacteria: DOM is lost!
   With bacteria: Energy recycles back into food web!
EOF
echo ""

echo -e "${MAGENTA}${BOLD}INTERPRETING YOUR MAG DATA:${NC}\n"

echo -e "${GREEN}📊 DIVERSITY PATTERNS:${NC}"
echo ""
echo "   HIGH RICHNESS (many species):"
echo "   ✓ Spatially heterogeneous environment"
echo "   ✓ Stable conditions (low disturbance)"
echo "   ✓ Resource diversity"
echo "   ✓ Example: Coastal, temperate, DCM"
echo ""
echo "   LOW RICHNESS (few species):"
echo "   ✓ Extreme conditions (hot, cold, acidic)"
echo "   ✓ Recent disturbance (mixing, bloom crash)"
echo "   ✓ Resource limitation"
echo "   ✓ Example: Bloom peak, hydrothermal vents"
echo ""

echo -e "${GREEN}📊 ABUNDANCE DISTRIBUTIONS:${NC}"
echo ""
echo "   EVEN DISTRIBUTION:"
echo "   → Many species, similar abundances"
echo "   → High competition, niche partitioning"
echo "   → Stable, mature community"
echo ""
echo "   SKEWED DISTRIBUTION:"
echo "   → Few dominant, many rare"
echo "   → Recent disturbance or bloom"
echo "   → One species 'winning'"
echo ""

echo -e "${GREEN}📊 FUNCTIONAL REDUNDANCY:${NC}"
echo ""
echo "   Multiple species, same function:"
echo "   ✓ Ecosystem resilience"
echo "   ✓ Insurance effect"
echo "   ✓ Different niches, same role"
echo ""
echo "   Example: Nitrogen fixation"
echo "   ├─ Trichodesmium (day)"
echo "   ├─ Crocosphaera (night)"
echo "   └─ Same function, different timing!"
echo ""

echo -e "${MAGENTA}${BOLD}METABOLIC GUILDS IN YOUR MAGS:${NC}\n"

cat << "EOF"
╔════════════════════════════════════════════════════════════╗
║                                                            ║
║  🌞 PRIMARY PRODUCERS                                      ║
║     └─ Cyanobacteria (photosynthesis: CO2 → organic C)   ║
║                                                            ║
║  ♻️  HETEROTROPHS                                          ║
║     ├─ SAR11 (organic carbon consumers)                   ║
║     ├─ Flavobacteria (particle degraders)                 ║
║     └─ Bacteroidetes (polymer specialists)                ║
║                                                            ║
║  ⚗️  CHEMOAUTOTROPHS                                       ║
║     ├─ Nitrifiers (ammonia → nitrite → nitrate)          ║
║     ├─ Sulfur oxidizers (H2S → SO4)                       ║
║     └─ Methane oxidizers (CH4 → CO2)                      ║
║                                                            ║
║  🔁 NITROGEN CYCLERS                                       ║
║     ├─ Fixers (N2 → NH3) - Trichodesmium                 ║
║     ├─ Nitrifiers (NH3 → NO3) - Archaea, Nitrospirae     ║
║     ├─ Denitrifiers (NO3 → N2) - Pseudomonas             ║
║     └─ Anammox (NH4 + NO2 → N2) - Planctomycetes         ║
║                                                            ║
║  💎 SPECIALISTS                                            ║
║     ├─ Iron oxidizers (Fe2+ → Fe3+)                       ║
║     ├─ Sulfate reducers (SO4 → H2S)                       ║
║     └─ Methanogens (CO2 + H2 → CH4)                       ║
║                                                            ║
╚════════════════════════════════════════════════════════════╝
EOF

echo ""
echo -e "${MAGENTA}${BOLD}ECOLOGICAL QUESTIONS TO ASK:${NC}\n"

echo "🤔 Is community structure driven by environment or history?"
echo "   → Compare similar environments: Same community?"
echo ""
echo "🤔 Are there keystone species?"
echo "   → Nitrogen fixers in N-limited waters?"
echo "   → Sulfur oxidizers at hydrothermal vents?"
echo ""
echo "🤔 Is functional redundancy high?"
echo "   → Multiple MAGs with same pathway?"
echo "   → Insurance against perturbation!"
echo ""
echo "🤔 Are there tightly linked species?"
echo "   → Co-occurrence networks"
echo "   → Metabolic dependencies (syntrophy)"
echo ""
echo "🤔 How does diversity change with depth?"
echo "   → Surface: High productivity, lower diversity"
echo "   → DCM: Highest diversity"
echo "   → Deep: Lower diversity, specialists"
echo ""

echo -e "${MAGENTA}${BOLD}COMMON ECOLOGICAL PITFALLS:${NC}\n"

echo "❌ Confusing presence with activity"
echo "   → DNA = potential, RNA = activity!"
echo ""
echo "❌ Ignoring rare taxa"
echo "   → They're the seed bank for future blooms"
echo ""
echo "❌ Assuming dominance = importance"
echo "   → Rare nitrogen fixers can control ecosystem!"
echo ""
echo "❌ Forgetting about viruses"
echo "   → Viral lysis releases nutrients"
echo "   → Virome sequencing is worth it!"
echo ""
echo "❌ Overlooking spatial scale"
echo "   → Micrometers: Single cell"
echo "   → Millimeters: Particle microenvironments"
echo "   → Meters: Water masses"
echo "   → Kilometers: Ocean provinces"
echo ""

echo -e "${MAGENTA}${BOLD}ADVANCED CONCEPTS:${NC}\n"

echo -e "${GREEN}🧬 HORIZONTAL GENE TRANSFER${NC}"
echo "   Bacteria share genes like files!"
echo "   → Antibiotic resistance spreads"
echo "   → Metabolic innovations transfer"
echo "   → Viruses as gene shuttles"
echo ""

echo -e "${GREEN}⚡ METABOLIC HANDOFFS${NC}"
echo "   Species A produces → Species B consumes"
echo "   Example:"
echo "   ├─ Diatoms die → Release DOM"
echo "   ├─ Flavobacteria → Degrade DOM"
echo "   ├─ Release CO2 → Cyanobacteria use it"
echo "   └─ Circle of life!"
echo ""

echo -e "${GREEN}🎯 RESOURCE RATIO THEORY${NC}"
echo "   N:P ratio determines winners"
echo "   ├─ High N:P → P-limited → Diatoms"
echo "   ├─ Low N:P → N-limited → Nitrogen fixers"
echo "   └─ Check nutrient data!"
echo ""

echo -e "${MAGENTA}${BOLD}RECOMMENDED ANALYSES FOR YOUR MAGS:${NC}\n"

echo "1. 🧮 Calculate diversity indices"
echo "   Shannon, Simpson, Pielou evenness"
echo "   → In R: vegan package"
echo ""
echo "2. 📊 Ordination (PCA, NMDS)"
echo "   Plot samples by community composition"
echo "   → Cluster by environment, not just geography!"
echo ""
echo "3. 🕸️  Co-occurrence networks"
echo "   Which MAGs co-occur?"
echo "   → Potential interactions, metabolic dependencies"
echo ""
echo "4. 🧬 Functional annotation"
echo "   What pathways are present?"
echo "   → KEGG, COG, Pfam databases"
echo ""
echo "5. 📈 Rarefaction curves"
echo "   Did you sequence enough?"
echo "   → Plateaus = good coverage"
echo ""

cat << "EOF"

     ╔═══════════════════════════════════════════════════╗
     ║  🦠 "Microbes are small, but their impact        ║
     ║      is planetary. Study them ecologically:      ║
     ║      Who's there, what are they doing,           ║
     ║      and why does it matter?"                    ║
     ╚═══════════════════════════════════════════════════╝

           Dr. Petra Microbe, Microbial Ecologist

EOF

echo -e "${GREEN}${BOLD}Need ecology insights? Ask me about:${NC}"
echo "• Community diversity & structure"
echo "• Metabolic guilds & interactions"
echo "• Ecological theory applications"
echo "• Statistical analyses"
echo ""
echo -e "${MAGENTA}🦠 Think ecologically! 🦠${NC}\n"
