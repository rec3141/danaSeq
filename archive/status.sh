#!/usr/bin/env bash
################################################################################
#                                                                              #
#  📊  DANA PIPELINE STATUS CHECKER  📊                                        #
#                                                                              #
#  Quick overview of what's installed and ready to roll!                      #
#                                                                              #
################################################################################

# Colors
readonly GREEN='\033[0;32m'
readonly RED='\033[0;31m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly BOLD='\033[1m'
readonly NC='\033[0m'

# Check if a command exists
cmd_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Print status
print_status() {
    local name="$1"
    local cmd="$2"
    local width=30

    printf "  %-${width}s " "$name"

    if cmd_exists "$cmd"; then
        local version=$($cmd --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1 || echo "")
        echo -e "${GREEN}✓${NC} ${version}"
    else
        echo -e "${RED}✗${NC} Not found"
    fi
}

echo ""
echo -e "${BLUE}${BOLD}╔═══════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}${BOLD}║                                                       ║${NC}"
echo -e "${BLUE}${BOLD}║        🔍 DANA PIPELINE STATUS CHECK 🔍               ║${NC}"
echo -e "${BLUE}${BOLD}║                                                       ║${NC}"
echo -e "${BLUE}${BOLD}╚═══════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${YELLOW}${BOLD}📦 CORE DEPENDENCIES${NC}"
echo ""

print_status "Bash" "bash"
print_status "Python 3" "python3"
print_status "R" "R"
print_status "AWK" "awk"
print_status "Parallel (GNU)" "parallel"

echo ""
echo -e "${YELLOW}${BOLD}🧬 SEQUENCING & ASSEMBLY${NC}"
echo ""

print_status "Flye assembler" "flye"
print_status "minimap2" "minimap2"
print_status "samtools" "samtools"
print_status "Racon" "racon"
print_status "Medaka" "medaka"

echo ""
echo -e "${YELLOW}${BOLD}🏷️ TAXONOMIC CLASSIFICATION${NC}"
echo ""

print_status "Kraken2" "kraken2"
print_status "Kaiju" "kaiju"
print_status "Prokka" "prokka"

echo ""
echo -e "${YELLOW}${BOLD}🗂️ BINNING TOOLS${NC}"
echo ""

print_status "SemiBin2" "semibin2"
print_status "MetaBAT2" "metabat2"
print_status "MaxBin2" "run_MaxBin.pl"
print_status "DAS Tool" "DAS_Tool"
print_status "CheckM2" "checkm2"

echo ""
echo -e "${YELLOW}${BOLD}🧹 PREPROCESSING${NC}"
echo ""

print_status "BBTools (bbduk)" "bbduk.sh"
print_status "Filtlong" "filtlong"
print_status "seqtk" "seqtk"

echo ""
echo -e "${YELLOW}${BOLD}💾 DATA MANAGEMENT${NC}"
echo ""

if cmd_exists "R"; then
    echo -n "  DuckDB (R package)             "
    if R --slave -e 'library(duckdb)' >/dev/null 2>&1; then
        echo -e "${GREEN}✓${NC}"
    else
        echo -e "${RED}✗${NC} Not installed"
    fi

    echo -n "  tidyverse (R package)          "
    if R --slave -e 'library(tidyverse)' >/dev/null 2>&1; then
        echo -e "${GREEN}✓${NC}"
    else
        echo -e "${RED}✗${NC} Not installed"
    fi

    echo -n "  leaflet (R package)            "
    if R --slave -e 'library(leaflet)' >/dev/null 2>&1; then
        echo -e "${GREEN}✓${NC}"
    else
        echo -e "${RED}✗${NC} Not installed"
    fi
else
    echo "  R not found - skipping R packages"
fi

echo ""
echo -e "${YELLOW}${BOLD}📁 DIRECTORY STRUCTURE${NC}"
echo ""

check_dir() {
    local dir="$1"
    local desc="$2"
    printf "  %-30s " "$desc"
    if [[ -d "$dir" ]]; then
        local count=$(find "$dir" -maxdepth 1 -type f -name "*.sh" -o -name "*.r" -o -name "*.R" -o -name "*.py" 2>/dev/null | wc -l | tr -d ' ')
        echo -e "${GREEN}✓${NC} ($count scripts)"
    else
        echo -e "${RED}✗${NC} Missing"
    fi
}

check_dir "nanopore_live" "Real-time processing"
check_dir "nanopore_mag" "MAG assembly"
check_dir "archive" "Archive"

echo ""
echo -e "${YELLOW}${BOLD}🎯 RECOMMENDED ACTIONS${NC}"
echo ""

# Check for missing critical tools
MISSING_CRITICAL=0

if ! cmd_exists "flye"; then
    echo -e "  ${RED}•${NC} Install Flye: ${BLUE}conda install -c bioconda flye${NC}"
    MISSING_CRITICAL=1
fi

if ! cmd_exists "kraken2"; then
    echo -e "  ${RED}•${NC} Install Kraken2: ${BLUE}conda install -c bioconda kraken2${NC}"
    MISSING_CRITICAL=1
fi

if ! cmd_exists "minimap2"; then
    echo -e "  ${RED}•${NC} Install minimap2: ${BLUE}conda install -c bioconda minimap2${NC}"
    MISSING_CRITICAL=1
fi

if ! cmd_exists "semibin2"; then
    echo -e "  ${RED}•${NC} Install SemiBin2: ${BLUE}pip install semibin${NC}"
    MISSING_CRITICAL=1
fi

if ! cmd_exists "metabat2"; then
    echo -e "  ${RED}•${NC} Install MetaBAT2: ${BLUE}conda install -c bioconda metabat2${NC}"
    MISSING_CRITICAL=1
fi

if [[ $MISSING_CRITICAL -eq 0 ]]; then
    echo -e "  ${GREEN}✓${NC} All critical tools installed! You're ready to sequence! 🎉"
fi

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  🌊 Ready to decode some oceans? 🌊              ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════╝${NC}"
echo ""
