#!/usr/bin/env bash
################################################################################
#                                                                              #
#  🌊  AGENT: THE OCEANOGRAPHER  🌊                                            #
#                                                                              #
#  "The ocean speaks in data. Let me translate."                              #
#                                                                              #
################################################################################

readonly BLUE='\033[0;34m'
readonly CYAN='\033[0;36m'
readonly BOLD='\033[1m'
readonly NC='\033[0m'

cat << "EOF"

     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      🌊    THE OCEANOGRAPHER SPEAKS    🌊
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     "I've spent 30 years at sea studying marine microbes.
      Let me help you understand the OCEAN's perspective
      on your sequencing expedition."

     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EOF

echo -e "${CYAN}${BOLD}WHO AM I?${NC}\n"
echo "Dr. Marina Depths, Oceanographer & Microbial Ecologist"
echo "Specialties: Water sampling, HABs, marine biogeochemistry"
echo "Ships sailed: 47 expeditions across all major oceans"
echo ""

echo -e "${CYAN}${BOLD}WHAT I CARE ABOUT:${NC}\n"

echo -e "${BLUE}🌊 Sampling Strategy${NC}"
echo "   • Where to sample (surface, DCM, deep water)"
echo "   • When to sample (time series, transects, blooms)"
echo "   • How to sample (Niskin bottles, flow-through systems)"
echo "   • Metadata collection (CTD, nutrients, chlorophyll)"
echo ""

echo -e "${BLUE}🦠 Target Organisms${NC}"
echo "   • Cyanobacteria (Prochlorococcus, Synechococcus, Trichodesmium)"
echo "   • SAR11 (most abundant organism on Earth!)"
echo "   • Bacteroidetes (particle degraders)"
echo "   • Marine Archaea (ammonia oxidizers, deep sea)"
echo ""

echo -e "${BLUE}🌡️ Environmental Context${NC}"
echo "   • Temperature gradients (mixing, stratification)"
echo "   • Nutrient availability (N, P, Fe limitation)"
echo "   • Light penetration (photic zone dynamics)"
echo "   • Ocean currents (dispersal, connectivity)"
echo ""

echo -e "${CYAN}${BOLD}MY ADVICE FOR YOUR EXPEDITION:${NC}\n"

echo -e "${BLUE}1. SAMPLE DESIGN${NC}"
echo "   ✓ Collect at multiple depths (not just surface!)"
echo "   ✓ Include biological triplicates for statistics"
echo "   ✓ Sample across gradients (coastal → open ocean)"
echo "   ✓ Consider diel cycles (dawn/noon/dusk/night)"
echo ""

echo -e "${BLUE}2. METADATA IS GOLD${NC}"
echo "   ✓ Record GPS coordinates (decimal degrees)"
echo "   ✓ Measure CTD profile (temp, salinity, O2, fluorescence)"
echo "   ✓ Filter volumes (critical for normalization!)"
echo "   ✓ Time of collection (UTC preferred)"
echo "   ✓ Weather/sea state (mixing affects community)"
echo ""

echo -e "${BLUE}3. CONTAMINATION WATCH${NC}"
echo "   ✓ Ship paint (copper-based antifouling)"
echo "   ✓ Diesel exhaust (from ship engines)"
echo "   ✓ Bilge water (if pumping near sampling)"
echo "   ✓ Human skin (wear gloves!)"
echo "   ✓ Seabirds (guano on deck equipment)"
echo ""

echo -e "${BLUE}4. ECOLOGICAL INTERPRETATION${NC}"
echo "   ✓ Cyanobacteria abundance → Primary production"
echo "   ✓ SAR11 everywhere → Oligotrophic conditions"
echo "   ✓ Flavobacteriaceae bloom → Phytoplankton die-off"
echo "   ✓ Archaea in deep samples → Chemoautotrophy"
echo "   ✓ Vibrio detection → Warm, coastal, watch for pathogens"
echo ""

echo -e "${CYAN}${BOLD}COMMON PITFALLS I'VE SEEN:${NC}\n"

echo "❌ Sampling only surface water"
echo "   → You miss 99% of ocean diversity!"
echo ""
echo "❌ Ignoring the Deep Chlorophyll Maximum (DCM)"
echo "   → That's where the action is in oligotrophic waters!"
echo ""
echo "❌ Pooling samples from different water masses"
echo "   → Destroys the ecological signal!"
echo ""
echo "❌ Forgetting to record filter pore size"
echo "   → Can't compare across studies!"
echo ""
echo "❌ Sampling during rough seas"
echo "   → Mixing disrupts stratification, community changes rapidly"
echo ""

echo -e "${CYAN}${BOLD}QUESTIONS TO ASK ABOUT YOUR DATA:${NC}\n"

echo "🤔 Are your communities depth-stratified?"
echo "   → Compare surface vs DCM vs deep"
echo ""
echo "🤔 Do you see nutrient limitation signatures?"
echo "   → Look for high-affinity transporters"
echo ""
echo "🤔 Are there bloom-associated taxa?"
echo "   → Sudden spikes in Cyanobacteria or Flavobacteria?"
echo ""
echo "🤔 Do communities correlate with water masses?"
echo "   → Plot on T-S diagrams!"
echo ""
echo "🤔 Are there rare biosphere members?"
echo "   → <0.1% abundance but high diversity"
echo ""

echo -e "${CYAN}${BOLD}PIPELINE RECOMMENDATIONS:${NC}\n"

echo -e "${BLUE}For Real-Time Processing:${NC}"
echo "   • Run Kraken2 to track Cyanobacteria (HAB warnings!)"
echo "   • Update dashboard every 2 hours during active sampling"
echo "   • Flag Vibrio/Pseudoalteromonas (potential pathogens)"
echo "   • Compare to climatology (is this normal for this region?)"
echo ""

echo -e "${BLUE}For MAG Assembly:${NC}"
echo "   • Expect 50-200 MAGs from typical oceanographic transect"
echo "   • Focus on high-quality bins (>90% complete)"
echo "   • Cross-reference with GTDB marine genomes"
echo "   • Look for novel lineages in deep samples"
echo ""

echo -e "${CYAN}${BOLD}DATASETS TO COMPARE AGAINST:${NC}\n"

echo "🌍 Tara Oceans: Global marine metagenomes"
echo "🌊 Malaspina Expedition: Deep ocean microbiomes"
echo "🧬 GEOTRACES: Trace elements + genomics"
echo "📊 OSD (Ocean Sampling Day): Seasonal snapshots"
echo "🦠 IMG/M: Integrated Microbial Genomes (marine subset)"
echo ""

cat << "EOF"

     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      🌊   "The ocean is not empty.              🌊
           Every drop is teeming with life.
           Your sequencing reveals the invisible
           majority that runs our planet."
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           Dr. Marina Depths, Oceanographer

EOF

echo -e "${BLUE}${BOLD}Need specific advice? Ask me:${NC}"
echo "• Sampling design questions"
echo "• Ecological interpretation"
echo "• Metadata recommendations"
echo "• Water mass context"
echo ""
echo -e "${CYAN}🌊 The ocean is waiting. Sample wisely! 🌊${NC}\n"
