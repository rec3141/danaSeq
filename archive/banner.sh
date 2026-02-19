#!/usr/bin/env bash
################################################################################
#                                                                              #
#  🎨  DANA PIPELINE BANNER  🎨                                                #
#                                                                              #
#  Show a beautiful welcome message for the pipeline!                          #
#                                                                              #
################################################################################

# Colors
readonly CYAN='\033[0;36m'
readonly BLUE='\033[0;34m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly MAGENTA='\033[0;35m'
readonly RED='\033[0;31m'
readonly BOLD='\033[1m'
readonly NC='\033[0m'

# Get terminal width for centering
TERM_WIDTH=$(tput cols 2>/dev/null || echo 80)

# Center text function
center_text() {
    local text="$1"
    local text_length=${#text}
    local padding=$(( (TERM_WIDTH - text_length) / 2 ))
    printf "%${padding}s" ""
    echo "$text"
}

clear

echo -e "${CYAN}"
cat << "EOF"
╔═══════════════════════════════════════════════════════════════════════════════╗
║                                                                               ║
║   ███╗   ███╗███████╗████████╗ █████╗  ██████╗ ███████╗███╗   ██╗ ██████╗   ║
║   ████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝████╗  ██║██╔═══██╗  ║
║   ██╔████╔██║█████╗     ██║   ███████║██║  ███╗█████╗  ██╔██╗ ██║██║   ██║  ║
║   ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██╔══╝  ██║╚██╗██║██║   ██║  ║
║   ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝███████╗██║ ╚████║╚██████╔╝  ║
║   ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝   ║
║                                                                               ║
║              🌊 OXFORD NANOPORE EDNA ANALYSIS PIPELINE 🌊                    ║
║                                                                               ║
║          Real-time metagenomic sequencing for oceanographic research         ║
║                                                                               ║
╚═══════════════════════════════════════════════════════════════════════════════╝
EOF
echo -e "${NC}"

echo -e "${BLUE}"
cat << "EOF"
                     🦠      🧬      🦠      🧬      🦠
                  ╔═══════════════════════════════════════╗
                  ║  ATCGATCGATCGATCGATCGATCGATCGATCGATCG ║
                  ║  TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG ║
                  ╚═══════════════════════════════════════╝
                     🔬   DECODE THE OCEANS   🔬
EOF
echo -e "${NC}\n"

echo -e "${GREEN}${BOLD}📂 PIPELINE STRUCTURE${NC}\n"

echo -e "${CYAN}  nanopore_live/${NC}     ⚡ Live analysis at sea"
echo -e "    ├─ Process reads as they stream"
echo -e "    ├─ Kraken2 taxonomic classification"
echo -e "    ├─ Prokka gene annotation"
echo -e "    └─ DuckDB integration & visualization"
echo ""

echo -e "${MAGENTA}  nanopore_mag/${NC}           🧬 Genome reconstruction"
echo -e "    ├─ Co-assembly with Flye"
echo -e "    ├─ Multi-tool binning (3-way consensus)"
echo -e "    ├─ Polishing with Racon + Medaka"
echo -e "    └─ Quality assessment & taxonomy"
echo ""

echo -e "${YELLOW}  archive/${NC}                💀 Deprecated scripts"
echo ""

echo -e "${GREEN}${BOLD}🚀 QUICK START${NC}\n"

echo -e "${CYAN}  Real-time Processing:${NC}"
echo -e "    cd nanopore_live"
echo -e "    ./24_process_reads_optimized.sh -i <input_dir> -K -P -S"
echo ""

echo -e "${MAGENTA}  MAG Assembly:${NC}"
echo -e "    cd nanopore_mag"
echo -e "    ./61_map_and_bin_optimized.sh"
echo ""

echo -e "${GREEN}${BOLD}🎯 KEY FEATURES${NC}\n"

echo -e "  ⚡ ${CYAN}Real-time${NC}       Process as sequencer streams data"
echo -e "  🎯 ${CYAN}Multi-tool${NC}      3-way binning consensus for robustness"
echo -e "  ✅ ${CYAN}Quality${NC}         QC at every step of the pipeline"
echo -e "  🚢 ${CYAN}Expedition${NC}      Optimized for shipboard deployment"
echo -e "  📊 ${CYAN}Visualization${NC}   Interactive maps & clustering"
echo -e "  🧬 ${CYAN}High-quality${NC}    Publication-ready MAGs"
echo -e "  🌊 ${CYAN}eDNA${NC}            Marine & freshwater specialized"
echo ""

echo -e "${GREEN}${BOLD}🎓 TARGET APPLICATIONS${NC}\n"

echo -e "  🌊 Marine microbial ecology"
echo -e "  🦠 Harmful algal bloom monitoring"
echo -e "  🚨 Waterborne pathogen surveillance"
echo -e "  🌍 Biodiversity assessments"
echo -e "  🧬 Environmental DNA profiling"
echo ""

echo -e "${GREEN}${BOLD}📚 DOCUMENTATION${NC}\n"

echo -e "  📖 Main README:               ${YELLOW}./README.md${NC}"
echo -e "  📖 Real-time processing:      ${YELLOW}./nanopore_live/README.md${NC}"
echo -e "  📖 MAG assembly:              ${YELLOW}./nanopore_mag/README.md${NC}"
echo ""

echo -e "${GREEN}${BOLD}🌍 ACTIVE EXPEDITIONS${NC}\n"

echo -e "  🚢 ${BLUE}CMO2025${NC} - California to Mexico Oceanographic Survey"
echo -e "  🧊 ${BLUE}QEI2025${NC} - Queen Elizabeth Islands Arctic Expedition"
echo ""

echo -e "${CYAN}"
cat << "EOF"
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║           🌊 DECODE THE OCEANS, ONE READ AT A TIME 🌊        ║
║                                                              ║
║                    🦠 → 🧬 → 💻 → 📊 → 🌍                    ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
EOF
echo -e "${NC}\n"

echo -e "${YELLOW}${BOLD}Now go forth and sequence! 🚀🔬🌊${NC}\n"
