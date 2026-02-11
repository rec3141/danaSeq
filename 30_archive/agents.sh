#!/usr/bin/env bash
################################################################################
#                                                                              #
#  🎭  DANA PIPELINE AGENTS  🎭                                                #
#                                                                              #
#  Your team of expert advisors!                                              #
#                                                                              #
################################################################################

readonly CYAN='\033[0;36m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly MAGENTA='\033[0;35m'
readonly BLUE='\033[0;34m'
readonly BOLD='\033[1m'
readonly NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AGENTS_DIR="${SCRIPT_DIR}/agents"

show_menu() {
    clear
    echo ""
    echo -e "${CYAN}╔═══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}║                                                               ║${NC}"
    echo -e "${CYAN}║           🎭  DANA PIPELINE EXPERT AGENTS  🎭                 ║${NC}"
    echo -e "${CYAN}║                                                               ║${NC}"
    echo -e "${CYAN}║           Your team of specialized advisors!                  ║${NC}"
    echo -e "${CYAN}║                                                               ║${NC}"
    echo -e "${CYAN}╚═══════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "${BOLD}Choose an expert to consult:${NC}"
    echo ""

    echo -e "${BLUE}  1) 🌊  The Oceanographer${NC}"
    echo "     Dr. Marina Depths - Sampling strategy, water masses, ecology"
    echo ""

    echo -e "${GREEN}  2) 💻  The Bioinformatician${NC}"
    echo "     Dr. Ada Pipeline - Code optimization, debugging, HPC"
    echo ""

    echo -e "${CYAN}  3) 🌊  The Ocean${NC}"
    echo "     The Ocean itself - Deep wisdom, planetary perspective"
    echo ""

    echo -e "${MAGENTA}  4) 🦠  The Microbial Ecologist${NC}"
    echo "     Dr. Petra Microbe - Community ecology, metabolic guilds"
    echo ""

    echo -e "${YELLOW}  5) 🎨  ALL AGENTS${NC}"
    echo "     Consult the entire team at once!"
    echo ""

    echo -e "${CYAN}  0) Exit${NC}"
    echo ""
    echo -n "Select (0-5): "
}

run_agent() {
    local agent_file="$1"

    if [[ ! -f "${agent_file}" ]]; then
        echo -e "${RED}Error: Agent script not found: ${agent_file}${NC}"
        return 1
    fi

    chmod +x "${agent_file}"
    bash "${agent_file}"

    echo ""
    echo -e "${CYAN}Press Enter to return to menu...${NC}"
    read
}

run_all_agents() {
    local agents=(
        "${AGENTS_DIR}/oceanographer.sh"
        "${AGENTS_DIR}/bioinformatician.sh"
        "${AGENTS_DIR}/ocean.sh"
        "${AGENTS_DIR}/microbial_ecologist.sh"
    )

    echo -e "${YELLOW}${BOLD}CONSULTING ALL AGENTS...${NC}\n"
    sleep 1

    for agent in "${agents[@]}"; do
        if [[ -f "${agent}" ]]; then
            chmod +x "${agent}"
            bash "${agent}"
            echo ""
            echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
            echo ""
            sleep 2
        fi
    done

    echo ""
    echo -e "${GREEN}${BOLD}✨ ALL AGENTS HAVE SPOKEN! ✨${NC}"
    echo ""
    echo -e "${CYAN}Press Enter to return to menu...${NC}"
    read
}

main() {
    # Make sure agents directory exists
    if [[ ! -d "${AGENTS_DIR}" ]]; then
        echo -e "${RED}Error: Agents directory not found: ${AGENTS_DIR}${NC}"
        exit 1
    fi

    while true; do
        show_menu
        read -r choice

        case $choice in
            1)
                run_agent "${AGENTS_DIR}/oceanographer.sh"
                ;;
            2)
                run_agent "${AGENTS_DIR}/bioinformatician.sh"
                ;;
            3)
                run_agent "${AGENTS_DIR}/ocean.sh"
                ;;
            4)
                run_agent "${AGENTS_DIR}/microbial_ecologist.sh"
                ;;
            5)
                run_all_agents
                ;;
            0)
                echo ""
                echo -e "${GREEN}Thanks for consulting the agents! 🎭${NC}"
                echo -e "${CYAN}May your sequencing be fruitful! 🧬🌊${NC}"
                echo ""
                exit 0
                ;;
            *)
                echo ""
                echo -e "${YELLOW}Invalid choice. Please select 0-5.${NC}"
                sleep 2
                ;;
        esac
    done
}

# Check if running directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
