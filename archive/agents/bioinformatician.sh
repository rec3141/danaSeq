#!/usr/bin/env bash
################################################################################
#                                                                              #
#  💻  AGENT: THE BIOINFORMATICIAN  💻                                         #
#                                                                              #
#  "Pipelines are poetry. Let me optimize yours."                             #
#                                                                              #
################################################################################

readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BOLD='\033[1m'
readonly NC='\033[0m'

cat << "EOF"

     ┌────────────────────────────────────────────────┐
     │  💻   THE BIOINFORMATICIAN SPEAKS   💻        │
     └────────────────────────────────────────────────┘

     "I've wrangled a petabyte of sequencing data and
      lived to tell the tale. Let me save you from my
      mistakes."

     ┌────────────────────────────────────────────────┐

EOF

echo -e "${GREEN}${BOLD}WHO AM I?${NC}\n"
echo "Dr. Ada Pipeline, Computational Biologist"
echo "Specialties: Metagenomics, HPC, workflow optimization"
echo "Lines of code written: Too many (most are Snakemake now)"
echo ""

echo -e "${GREEN}${BOLD}WHAT I CARE ABOUT:${NC}\n"

echo -e "${YELLOW}⚙️  Pipeline Efficiency${NC}"
echo "   • Minimize I/O bottlenecks"
echo "   • Parallelize everything possible"
echo "   • Cache intermediate results"
echo "   • Fail fast with clear errors"
echo ""

echo -e "${YELLOW}📊 Data Quality${NC}"
echo "   • GIGO: Garbage In, Garbage Out"
echo "   • QC at every step"
echo "   • Track contamination"
echo "   • Validate against known controls"
echo ""

echo -e "${YELLOW}🔄 Reproducibility${NC}"
echo "   • Version all software"
echo "   • Document parameters"
echo "   • Containerize workflows"
echo "   • Provide example datasets"
echo ""

echo -e "${GREEN}${BOLD}MY ADVICE FOR DANA PIPELINE:${NC}\n"

echo -e "${YELLOW}1. BEFORE YOU START${NC}"
echo "   ✓ Check dependencies: ./status.sh"
echo "   ✓ Test with small dataset first (100 reads)"
echo "   ✓ Estimate resources (RAM, disk, time)"
echo "   ✓ Set up logging directory"
echo "   ✓ Plan checkpointing strategy"
echo ""

echo -e "${YELLOW}2. REAL-TIME PROCESSING TIPS${NC}"
echo "   ✓ Use -K flag for Kraken (taxonomy is fast!)"
echo "   ✓ Skip Prokka until post-expedition (it's slow)"
echo "   ✓ Set --max-duration to avoid overruns"
echo "   ✓ Monitor DuckDB size (it grows fast!)"
echo "   ✓ Run dashboard in tmux/screen"
echo ""

echo -e "${YELLOW}3. MAG ASSEMBLY OPTIMIZATION${NC}"
echo "   ✓ Flye: Use --meta mode, --threads to max"
echo "   ✓ minimap2: -x map-ont for Nanopore"
echo "   ✓ Binning: Run 3 tools overnight"
echo "   ✓ CheckM2: MUCH faster than CheckM v1"
echo "   ✓ Polish only high-quality bins (>70% complete)"
echo ""

echo -e "${YELLOW}4. COMPUTE RESOURCES${NC}"
echo "   📦 Real-time processing:"
echo "      • RAM: 32GB minimum, 64GB comfortable"
echo "      • CPU: 16+ cores (Kraken loves threads)"
echo "      • Disk: 500GB for typical expedition"
echo "      • Time: ~2-4 hours per barcode"
echo ""
echo "   🧬 MAG assembly:"
echo "      • RAM: 128GB minimum, 256GB better"
echo "      • CPU: 32+ cores (Flye, MetaBAT parallel)"
echo "      • Disk: 1-2TB (assemblies are BIG)"
echo "      • Time: 1-3 days for complete pipeline"
echo ""

echo -e "${GREEN}${BOLD}COMMON ERRORS I'VE DEBUGGED 1000 TIMES:${NC}\n"

echo "🐛 \"Out of memory killed\""
echo "   Fix: Reduce threads, increase RAM, or split input"
echo ""
echo "🐛 \"Kraken database not found\""
echo "   Fix: Set KRAKEN_DB environment variable"
echo ""
echo "🐛 \"No such file or directory: assembly.fasta\""
echo "   Fix: Check previous step logs, assembly may have failed"
echo ""
echo "🐛 \"DAS_Tool: No bins passed quality filter\""
echo "   Fix: Lower --score_threshold, check bin quality"
echo ""
echo "🐛 \"Medaka: CUDA out of memory\""
echo "   Fix: Use CPU version or reduce batch size"
echo ""

echo -e "${GREEN}${BOLD}PERFORMANCE TUNING:${NC}\n"

echo -e "${YELLOW}❌ DON'T DO THIS:${NC}"
echo "   • Process all barcodes sequentially (use parallel!)"
echo "   • Re-run entire pipeline after failure (use checkpoints!)"
echo "   • Store everything in one directory (organize!)"
echo "   • Run on login node (use compute nodes!)"
echo "   • Ignore errors (check logs!)"
echo ""

echo -e "${YELLOW}✅ DO THIS:${NC}"
echo "   • Use GNU parallel for barcode processing"
echo "   • Check for existing outputs before rerunning"
echo "   • Separate by project/sample/date"
echo "   • Submit to SLURM/PBS queue"
echo "   • grep ERROR logs/ after each run"
echo ""

echo -e "${GREEN}${BOLD}QUALITY CONTROL CHECKLIST:${NC}\n"

echo "📊 After Real-Time Processing:"
echo "   □ Read count per sample >10K"
echo "   □ Mean read length >2kb"
echo "   □ Mean read quality >Q10"
echo "   □ Kraken classified >50%"
echo "   □ No single taxon >90% (likely contamination)"
echo "   □ DuckDB queryable"
echo ""

echo "🧬 After MAG Assembly:"
echo "   □ Assembly N50 >10kb"
echo "   □ Total contigs <100K"
echo "   □ Recovered >20 MAGs"
echo "   □ High-quality bins >5"
echo "   □ Contamination <10% for all bins"
echo "   □ Taxonomic assignment >80% of bins"
echo ""

echo -e "${GREEN}${BOLD}DEBUGGING WORKFLOW:${NC}\n"

echo "1. Check the logs!"
echo "   tail -100 logs/barcode01.log"
echo ""
echo "2. Verify inputs exist"
echo "   ls -lh input/barcode01/*.fastq.gz"
echo ""
echo "3. Test tool individually"
echo "   kraken2 --threads 4 --db \$KRAKEN_DB test.fq"
echo ""
echo "4. Check disk space"
echo "   df -h ."
echo ""
echo "5. Validate output format"
echo "   head -20 output/assembly.fasta"
echo ""

echo -e "${GREEN}${BOLD}ADVANCED TIPS:${NC}\n"

echo -e "${YELLOW}🚀 Speed Hacks:${NC}"
echo "   • Use ramdisk (/dev/shm) for temp files"
echo "   • Compress with pigz (parallel gzip)"
echo "   • Use samtools view -@ for parallel BAM reading"
echo "   • Cache Kraken2 database in RAM (--memory-mapping)"
echo ""

echo -e "${YELLOW}💾 Storage Hacks:${NC}"
echo "   • Delete intermediate BAMs after depth calculation"
echo "   • Compress assemblies with gzip -9"
echo "   • Symlink raw reads instead of copying"
echo "   • Archive to tape after publication"
echo ""

echo -e "${YELLOW}🔧 Troubleshooting Hacks:${NC}"
echo "   • Run with -d (debug mode) to see commands"
echo "   • Use strace to debug file access issues"
echo "   • htop to monitor CPU/RAM in real-time"
echo "   • iostat -x 1 to find I/O bottlenecks"
echo ""

echo -e "${GREEN}${BOLD}RECOMMENDED READING:${NC}\n"

echo "📚 Metagenomic assembly:"
echo "   • Nurk et al. 2017 (metaSPAdes)"
echo "   • Kolmogorov et al. 2020 (metaFlye)"
echo ""
echo "📚 Binning methods:"
echo "   • Pan et al. 2023 (SemiBin2)"
echo "   • Kang et al. 2019 (MetaBAT2)"
echo "   • Sieber et al. 2018 (DAS_Tool)"
echo ""
echo "📚 Best practices:"
echo "   • Bowers et al. 2017 (MIMAG standards)"
echo "   • Nayfach et al. 2021 (CheckM2)"
echo ""

cat << "EOF"

     ┌────────────────────────────────────────────────┐
     │  💻  "Bioinformatics is 90% data wrangling, │
     │       9% debugging, and 1% actual science.    │
     │       But that 1% changes the world."         │
     └────────────────────────────────────────────────┘

           Dr. Ada Pipeline, Bioinformatician

EOF

echo -e "${YELLOW}${BOLD}Need help? Ask me about:${NC}"
echo "• Pipeline optimization"
echo "• Resource requirements"
echo "• Error debugging"
echo "• Best practices"
echo ""
echo -e "${GREEN}💻 Code smarter, not harder! 💻${NC}\n"
