/*
 * fastq_filter — FASTQ dedup + quality/length filter.
 *
 * Replaces: cat *.fastq.gz | dedup | filtlong -t TARGET
 *
 * Pipe-friendly. With --target_bases the default is a two-pass bucket sort:
 * reads are spilled to score buckets in one decompression pass, then buckets
 * are emitted from the top down until the budget is met, so the kept subset
 * depends only on read scores and never on input order. --onepass restores the
 * legacy single-pass behaviour, which scores and immediately accepts or rejects
 * against a dynamically adjusted threshold converging on target_bases (biased
 * to undershoot for memory safety) and is therefore order-dependent.
 *
 * Scoring follows filtlong's algorithm (without global z-score normalization):
 *   length_score  = 100 * len / (len + 5000)
 *   mean_quality  = mean(1 - 10^(-Q/10)) * 100
 *   window_quality = min sliding-window avg * 100
 *   final = sqrt(length_score * mean_quality) * (0.5 + 0.5 * window_ratio)
 *
 * Compile: g++ -O2 -o fastq_filter fastq_filter.cpp -lz -lpthread
 * Usage:   fastq_filter -t 40000000000 input1.fastq.gz input2.fastq.gz > out.fastq
 *          cat *.fastq.gz | fastq_filter -t 40000000000 > out.fastq
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cerrno>
#include <memory>
#include <list>
#include <stdexcept>
#include <csignal>
#include <getopt.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/stat.h>
#include <sys/resource.h>

// Precomputed Phred+33 -> P(correct) lookup
static double PHRED_LUT[256];

static void io_check(bool ok, const char* operation) {
    if (!ok) throw std::runtime_error(std::string(operation) + ": " + strerror(errno));
}

// Private directory prevents collisions between concurrent processes on shared
// scratch. A small LRU bounds descriptors and buffering independently of buckets.
class SpillBuckets {
    std::string directory;
    std::vector<FILE*> files;
    std::list<size_t> recent;
    std::vector<std::list<size_t>::iterator> positions;
    size_t opened = 0;
    size_t max_open = 512;
public:
    std::vector<std::string> paths;
    explicit SpillBuckets(int count, const char* parent)
        : files(count, nullptr), positions(count), paths(count) {
        if (!count) return;
        struct rlimit limit;
        io_check(getrlimit(RLIMIT_NOFILE, &limit) == 0, "get descriptor limit");
        // Preserve the original default's buffering when possible, but leave
        // descriptors for input/output/runtime and support low-limit systems.
        max_open = std::min<rlim_t>(512, limit.rlim_cur > 16 ? limit.rlim_cur - 16 : 1);
        std::string pattern = std::string(parent) + "/fastq_filter.XXXXXX";
        std::vector<char> name(pattern.begin(), pattern.end());
        name.push_back('\0');
        io_check(mkdtemp(name.data()) != nullptr, "create spill directory");
        directory = name.data();
    }
    void close_one(size_t b) {
        FILE* fp = files[b];
        files[b] = nullptr;
        recent.erase(positions[b]);
        --opened;
        io_check(fclose(fp) == 0, "flush/close spill file");
    }
    FILE* writer(size_t b) {
        if (!files[b]) {
            if (opened == max_open) {
                close_one(recent.back());
            }
            if (paths[b].empty()) paths[b] = directory + "/" + std::to_string(b);
            files[b] = fopen(paths[b].c_str(), "ab");
            io_check(files[b] != nullptr, "open spill file");
            ++opened;
            recent.push_front(b);
            positions[b] = recent.begin();
            setvbuf(files[b], nullptr, _IOFBF, 1 << 16);
        }
        recent.splice(recent.begin(), recent, positions[b]);
        return files[b];
    }
    void finish() {
        for (size_t b = 0; b < files.size(); ++b) if (files[b]) close_one(b);
    }
    ~SpillBuckets() {
        for (auto fp : files) if (fp) fclose(fp);
        for (const auto& path : paths) if (!path.empty()) unlink(path.c_str());
        if (!directory.empty()) rmdir(directory.c_str());
    }
};

static void init_phred_lut() {
    for (int i = 0; i < 256; i++) PHRED_LUT[i] = 0.0;
    for (int q = 0; q < 94; q++) {
        PHRED_LUT[q + 33] = 1.0 - pow(10.0, -q / 10.0);
    }
}

struct ReadRecord {
    std::string header;
    std::string seq;
    std::string qual;
};

static inline double score_read(const char* qual, int length, int window_size) {
    if (length == 0) return 0.0;

    double length_score = 100.0 * length / (length + 5000.0);

    // Mean quality
    double total_q = 0.0;
    for (int i = 0; i < length; i++) {
        total_q += PHRED_LUT[(unsigned char)qual[i]];
    }
    double mean_q = total_q / length * 100.0;
    if (mean_q <= 0) return 0.0;

    // Window quality: min sliding window average
    double window_q;
    if (length <= window_size) {
        window_q = mean_q;
    } else {
        double wsum = 0.0;
        for (int i = 0; i < window_size; i++) {
            wsum += PHRED_LUT[(unsigned char)qual[i]];
        }
        double min_wsum = wsum;
        for (int i = 1; i <= length - window_size; i++) {
            wsum += PHRED_LUT[(unsigned char)qual[i + window_size - 1]]
                  - PHRED_LUT[(unsigned char)qual[i - 1]];
            if (wsum < min_wsum) min_wsum = wsum;
        }
        window_q = min_wsum / window_size * 100.0;
    }

    double window_ratio = (window_q < mean_q) ? window_q / mean_q : 1.0;
    return sqrt(length_score * mean_q) * (0.5 + 0.5 * window_ratio);
}

// Inverse normal CDF approximation (Beasley-Springer-Moro)
static double inv_normal(double p) {
    if (p <= 0.0) return -4.0;
    if (p >= 1.0) return 4.0;
    if (p <= 0.5) {
        double t = sqrt(-2.0 * log(p));
        return t - (2.515517 + 0.802853*t + 0.010328*t*t) /
               (1.0 + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t);
    } else {
        double t = sqrt(-2.0 * log(1.0 - p));
        return -(t - (2.515517 + 0.802853*t + 0.010328*t*t) /
               (1.0 + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t));
    }
}

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] [input.fastq.gz ...]\n"
        "\n  Output thresholds (compatible with filtlong):\n"
        "  -t, --target_bases N   Keep only the best reads up to this many total bases\n"
        "  -p, --keep_percent F   Keep only this percentage of the best reads (1-100)\n"
        "  --min_length N         Minimum length threshold (default: 0)\n"
        "  --min_mean_q F         Minimum mean quality threshold (0-100, default: 0)\n"
        "  --min_window_q F       Minimum window quality threshold (0-100, default: 0)\n"
        "\n  Additional options:\n"
        "  --window_size N        Quality scoring window (default: 250)\n"
        "  -D, --no_dedupe        Skip deduplication (filtlong has no dedup)\n"
        "  -o, --output FILE      Output file (default: stdout, .gz for gzipped)\n"
        "  -h, --help             Show this help\n"
        "\n  Selection mode (--target_bases only):\n"
        "  --twopass              Default. Spill reads to score buckets, then emit\n"
        "                         the best --target_bases. Output does not depend on\n"
        "                         input order. Needs scratch space for one\n"
        "                         uncompressed copy of the input.\n"
        "  --onepass              Legacy streaming threshold. Lower disk use, but the\n"
        "                         selected subset depends on input order.\n"
        "  --spill_dir DIR        Bucket directory (default: $SLURM_TMPDIR, else /tmp)\n"
        "  --buckets N            Score buckets for two-pass mode (default: 512)\n"
        "\nReads from stdin if no files given.\n"
        "Deduplication is ON by default (unique to fastq_filter).\n"
        "\n--onepass decides each read on arrival against a threshold built only from\n"
        "the reads seen so far, so it is lenient early and strict late: two groups\n"
        "with identical length and quality distributions keep 1518 vs 501 reads\n"
        "purely by position in the stream. Use it only where reproducibility and\n"
        "even sampling across inputs do not matter.\n"
        "\n--keep_percent always streams; it has no two-pass implementation.\n", prog);
}

static int run_main(int argc, char** argv) {
    init_phred_lut();

    long long target_bases = 0;
    double keep_percent = 0.0;  // 0 = disabled
    int min_length = 0;
    double min_mean_q = 0.0;
    double min_window_q = 0.0;
    int window_size = 250;
    bool dedupe = true;
    const char* output_path = nullptr;
    bool onepass = false;
    bool twopass_explicit = false;
    const char* spill_dir = nullptr;
    int n_buckets = 512;

    static struct option long_opts[] = {
        {"target_bases",  required_argument, 0, 't'},
        {"keep_percent",  required_argument, 0, 'p'},
        {"min_length",    required_argument, 0, 'm'},
        {"min_mean_q",    required_argument, 0, 'Q'},
        {"min_window_q",  required_argument, 0, 'W'},
        {"window_size",   required_argument, 0, 'w'},
        {"no_dedupe",     no_argument,       0, 'D'},
        {"output",        required_argument, 0, 'o'},
        {"onepass",       no_argument,       0, 'O'},
        {"twopass",       no_argument,       0, 'T'},
        {"spill_dir",     required_argument, 0, 'S'},
        {"buckets",       required_argument, 0, 'B'},
        {"help",          no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "t:p:m:w:Do:h", long_opts, nullptr)) != -1) {
        switch (c) {
            case 't': target_bases = atoll(optarg); break;
            case 'p': keep_percent = atof(optarg); break;
            case 'm': min_length = atoi(optarg); break;
            case 'Q': min_mean_q = atof(optarg); break;
            case 'W': min_window_q = atof(optarg); break;
            case 'w': window_size = atoi(optarg); break;
            case 'D': dedupe = false; break;
            case 'o': output_path = optarg; break;
            case 'O': onepass = true; break;
            case 'T': twopass_explicit = true; break;
            case 'S': spill_dir = optarg; break;
            case 'B': n_buckets = atoi(optarg); break;
            case 'h': usage(argv[0]); return 0;
            default:  usage(argv[0]); return 1;
        }
    }

    bool filter_by_bases = target_bases > 0;
    bool filter_by_pct = keep_percent > 0.0 && keep_percent < 100.0;
    bool filter_by_qual = min_mean_q > 0.0 || min_window_q > 0.0;
    bool filtering = filter_by_bases || filter_by_pct;

    if (window_size <= 0) {
        fprintf(stderr, "[fastq_filter] --window_size must be positive\n");
        return 1;
    }

    // Two-pass bucket selection is the DEFAULT wherever --target_bases applies,
    // because it is the only mode whose output does not depend on input order.
    // The single-pass streaming threshold remains available as --onepass: it
    // decides each read on arrival against a bar built from the reads seen so
    // far, so it is lenient early and strict late and the kept subset shifts
    // with input order. Percentage mode has no two-pass implementation and
    // always streams.
    bool twopass = filter_by_bases && !filter_by_pct && !onepass;

    if (onepass && twopass_explicit) {
        fprintf(stderr, "[fastq_filter] --onepass and --twopass are mutually exclusive\n");
        return 1;
    }
    if (twopass_explicit && !filter_by_bases) {
        fprintf(stderr, "[fastq_filter] --twopass requires --target_bases\n");
        return 1;
    }
    if (twopass_explicit && filter_by_pct) {
        fprintf(stderr, "[fastq_filter] --twopass does not support --keep_percent\n");
        return 1;
    }
    if (twopass) {
        if (n_buckets < 2 || n_buckets > 65536) {
            fprintf(stderr, "[fastq_filter] --buckets must be between 2 and 65536\n");
            return 1;
        }
        if (!spill_dir) {
            spill_dir = getenv("SLURM_TMPDIR");
            if (!spill_dir || !*spill_dir) spill_dir = "/tmp";
        }
    }

    // For percentage mode, compute the z-score threshold for the target percentile.
    // Assuming roughly normal score distribution:
    //   keep top K% → accept scores above the (100-K)th percentile
    //   threshold = mean + z * std, where z = inv_normal(1 - K/100)
    double pct_z = 0.0;
    if (filter_by_pct) {
        // Pre-compute z-score for percentage mode.
        // keep_percent=25 → keep top 25% → reject bottom 75%
        // → threshold at 75th percentile → z = inv_normal(0.75)
        pct_z = inv_normal(1.0 - keep_percent / 100.0);
    }

    // Collect input files
    std::vector<const char*> inputs;
    for (int i = optind; i < argc; i++) {
        inputs.push_back(argv[i]);
    }
    if (inputs.empty()) inputs.push_back("-");

    // Estimate total bases from file sizes
    long long est_total = 0;
    if (filtering) {
        for (auto& f : inputs) {
            if (strcmp(f, "-") == 0) { est_total = 0; break; }
            struct stat st;
            if (stat(f, &st) == 0) {
                // Conservative lower-bound estimate of seq bases.
                // gzip FASTQ: ~4x decompression, ~25% is sequence.
                // Use file_size * 0.8 as conservative floor.
                est_total += (long long)(st.st_size * 0.8);
            }
        }
    }

    // Output
    FILE* out = stdout;
    gzFile gz_out = nullptr;
    if (output_path) {
        size_t olen = strlen(output_path);
        if (olen > 3 && strcmp(output_path + olen - 3, ".gz") == 0) {
            gz_out = gzopen(output_path, "wb3");
            if (!gz_out) { fprintf(stderr, "Cannot open %s\n", output_path); return 1; }
        } else {
            out = fopen(output_path, "w");
            if (!out) { fprintf(stderr, "Cannot open %s\n", output_path); return 1; }
        }
    }

    // Dedup set
    std::unordered_set<std::string> seen_ids;

    // Welford online stats
    long long n_scored = 0;
    double score_mean = 0.0, score_m2 = 0.0;
    long long scored_bases = 0;
    double threshold = 0.0;

    // Stats
    long long total_reads = 0, dup_reads = 0, short_reads = 0, qual_reads = 0;
    long long accepted_reads = 0, accepted_bases = 0;

    // Score histograms (1000 bins from 0-100)
    static const int HIST_BINS = 1000;
    long long score_hist[HIST_BINS];      // all scored reads
    long long score_hist_kept[HIST_BINS]; // accepted reads
    memset(score_hist, 0, sizeof(score_hist));
    memset(score_hist_kept, 0, sizeof(score_hist_kept));

    // Line buffers
    const int MAXLINE = 1 << 20; // 1 MB max line
    char* hdr_buf = (char*)malloc(MAXLINE);
    char* seq_buf = (char*)malloc(MAXLINE);
    char* plus_buf = (char*)malloc(MAXLINE);
    char* qual_buf = (char*)malloc(MAXLINE);

    // --twopass spill buckets.
    //
    // Reads are written to one of n_buckets plain-FASTQ files chosen by score,
    // so the input is decompressed exactly once. Selection then walks buckets
    // from the top down until target_bases is reached. Which reads are kept
    // depends only on their scores, never on arrival order.
    //
    // Sort compact indices one bucket at a time; records remain on disk.
    SpillBuckets buckets(twopass ? n_buckets : 0, spill_dir);
    if (twopass)
        fprintf(stderr, "[fastq_filter] --twopass: %d buckets in %s\n",
                n_buckets, spill_dir);

    auto write_bytes = [&](const char* data, size_t length) {
        if (gz_out) io_check(gzwrite(gz_out, data, length) == (int)length, "write gzip output");
        else io_check(fwrite(data, 1, length, out) == length, "write output");
    };

    auto write_record = [&](const char* hdr, int hlen,
                            const char* seq, int slen,
                            const char* qual, int qlen) {
        write_bytes(hdr, hlen);
        write_bytes("\n", 1);
        write_bytes(seq, slen);
        write_bytes("\n+\n", 3);
        write_bytes(qual, qlen);
        write_bytes("\n", 1);
    };

    for (auto& fname : inputs) {
        gzFile gz;
        if (strcmp(fname, "-") == 0) {
            gz = gzdopen(fileno(stdin), "r");
        } else {
            gz = gzopen(fname, "r");
        }
        if (!gz) {
            throw std::runtime_error(std::string("Cannot open input ") + fname);
        }
        gzbuffer(gz, 1 << 18); // 256KB buffer

        while (true) {
            if (!gzgets(gz, hdr_buf, MAXLINE)) break;
            if (!gzgets(gz, seq_buf, MAXLINE) || !gzgets(gz, plus_buf, MAXLINE) ||
                !gzgets(gz, qual_buf, MAXLINE))
                throw std::runtime_error("Incomplete FASTQ record");
            if (!strchr(hdr_buf, '\n') || !strchr(seq_buf, '\n') ||
                !strchr(plus_buf, '\n') ||
                (!strchr(qual_buf, '\n') && !gzeof(gz)))
                throw std::runtime_error("FASTQ line exceeds buffer limit");

            total_reads++;

            // Strip newlines
            int hlen = strlen(hdr_buf);
            if (hlen > 0 && hdr_buf[hlen-1] == '\n') hdr_buf[--hlen] = '\0';
            int slen = strlen(seq_buf);
            if (slen > 0 && seq_buf[slen-1] == '\n') seq_buf[--slen] = '\0';
            int qlen = strlen(qual_buf);
            if (qlen > 0 && qual_buf[qlen-1] == '\n') qual_buf[--qlen] = '\0';
            if (hdr_buf[0] != '@' || plus_buf[0] != '+' || !slen || slen != qlen)
                throw std::runtime_error("Invalid FASTQ record or sequence/quality length mismatch");

            // Dedup by read ID
            if (dedupe && !twopass) {
                // Extract ID: skip '@', take up to first space/tab
                const char* start = hdr_buf + 1;
                const char* end = start;
                while (*end && *end != ' ' && *end != '\t') end++;
                std::string rid(start, end - start);
                if (!seen_ids.insert(rid).second) {
                    dup_reads++;
                    continue;
                }
            }

            // Min length
            if (slen < min_length) {
                short_reads++;
                continue;
            }

            // Min quality thresholds (filtlong-compatible, applied before scoring)
            if (filter_by_qual) {
                double total_q = 0.0;
                for (int i = 0; i < slen; i++)
                    total_q += PHRED_LUT[(unsigned char)qual_buf[i]];
                double mq = total_q / slen * 100.0;
                if (mq < min_mean_q) { qual_reads++; continue; }
                if (min_window_q > 0) {
                    double wq;
                    if (slen <= window_size) {
                        wq = mq;
                    } else {
                        double wsum = 0.0;
                        for (int i = 0; i < window_size; i++)
                            wsum += PHRED_LUT[(unsigned char)qual_buf[i]];
                        double min_wsum = wsum;
                        for (int i = 1; i <= slen - window_size; i++) {
                            wsum += PHRED_LUT[(unsigned char)qual_buf[i + window_size - 1]]
                                  - PHRED_LUT[(unsigned char)qual_buf[i - 1]];
                            if (wsum < min_wsum) min_wsum = wsum;
                        }
                        wq = min_wsum / window_size * 100.0;
                    }
                    if (wq < min_window_q) { qual_reads++; continue; }
                }
            }

            // No score-based filtering — just dedup/length/quality and output
            if (!filtering) {
                write_record(hdr_buf, hlen, seq_buf, slen, qual_buf, qlen);
                accepted_reads++;
                accepted_bases += slen;
                continue;
            }

            // Score
            double score = score_read(qual_buf, slen, window_size);

            // --twopass: spill to the score bucket and decide nothing yet.
            if (twopass) {
                int b = (int)(score / 100.0 * n_buckets);
                if (b >= n_buckets) b = n_buckets - 1;
                if (b < 0) b = 0;
                FILE* bf = buckets.writer(b);
                io_check(fwrite(hdr_buf, 1, hlen, bf) == (size_t)hlen &&
                         fputc('\n', bf) != EOF &&
                         fwrite(seq_buf, 1, slen, bf) == (size_t)slen &&
                         fputs("\n+\n", bf) != EOF &&
                         fwrite(qual_buf, 1, qlen, bf) == (size_t)qlen &&
                         fputc('\n', bf) != EOF, "write spill record");
                n_scored++;
                scored_bases += slen;
                continue;
            }

            // Welford online stats + histogram
            n_scored++;
            scored_bases += slen;
            int bin = (int)(score / 100.0 * HIST_BINS);
            if (bin >= HIST_BINS) bin = HIST_BINS - 1;
            if (bin < 0) bin = 0;
            score_hist[bin]++;
            double delta = score - score_mean;
            score_mean += delta / n_scored;
            score_m2 += delta * (score - score_mean);

            // Update threshold: every read for first 100 (fast convergence),
            // then every 100 reads (cheap maintenance)
            if (n_scored == 10 || (n_scored > 10 && n_scored % 100 == 0)) {
                double score_std = (n_scored > 1) ? sqrt(score_m2 / n_scored) : 0;
                if (score_std > 0) {
                    if (filter_by_pct) {
                        // Percentage mode: use histogram percentile with feedback.
                        // Track by read count (not bases) for accurate percentages.
                        double target_frac = keep_percent / 100.0;
                        double actual_frac = (n_scored > 0)
                            ? (double)accepted_reads / n_scored : 1.0;
                        // Error: positive = keeping too much
                        double error = actual_frac - target_frac;
                        // Compensate with gain=3: if we're 5pp over, aim 15pp stricter
                        double comp_frac = target_frac - 3.0 * error;
                        if (comp_frac < 0.01) comp_frac = 0.01;
                        if (comp_frac > 0.99) comp_frac = 0.99;

                        long long cutoff_count = (long long)(n_scored * (1.0 - comp_frac));
                        long long cumsum = 0;
                        threshold = 0.0;
                        for (int b = 0; b < HIST_BINS; b++) {
                            cumsum += score_hist[b];
                            if (cumsum >= cutoff_count) {
                                threshold = (b + 0.5) / HIST_BINS * 100.0;
                                break;
                            }
                        }
                    } else {
                        // Target bases mode
                        double remaining_need = (double)target_bases - accepted_bases;
                        double remaining_input = (double)scored_bases;
                        if (remaining_input < 1) remaining_input = 1;

                        // If target exceeds file-size estimate, accept all
                        if (est_total > 0 && target_bases >= est_total) {
                            threshold = 0.0;
                        } else if (remaining_need <= 0) {
                            threshold = score_mean + 2.0 * score_std;
                        } else {
                            double need_frac = remaining_need / remaining_input;
                            if (need_frac >= 1.0) {
                                threshold = 0.0;
                            } else {
                                if (need_frac < 0.01) need_frac = 0.01;
                                double z = inv_normal(1.0 - need_frac);
                                threshold = score_mean + z * score_std;
                            }
                        }
                    }
                }
            }

            // Accept/reject
            if (score >= threshold) {
                write_record(hdr_buf, hlen, seq_buf, slen, qual_buf, qlen);
                accepted_reads++;
                accepted_bases += slen;
                score_hist_kept[bin]++;
            }
        }

        int input_error = Z_OK;
        const char* input_message = gzerror(gz, &input_error);
        std::string message = input_message ? input_message : "input read failure";
        int close_error = gzclose(gz);
        if ((input_error != Z_OK && input_error != Z_STREAM_END) || close_error != Z_OK)
            throw std::runtime_error("Read input: " + message);
    }

    // Global greedy order: descending score, ID hash, then exact record bytes.
    // Every bucket is sorted, including fully accepted ones. Continue to lower
    // buckets when a read does not fit: bucket boundaries never affect selection.
    if (twopass) {
        buckets.finish(); // Detect buffered write failures before consuming spills.
        struct Cand {
            double score;
            uint64_t idhash;
            long off;
            int len;
        };
        for (int b = n_buckets - 1; b >= 0 && accepted_bases < target_bases; --b) {
            if (buckets.paths[b].empty()) continue;
            std::unique_ptr<FILE, decltype(&fclose)> file(
                fopen(buckets.paths[b].c_str(), "rb"), fclose);
            io_check(file != nullptr, "open spill for reading");
            FILE* bf = file.get();
            auto read_record = [&]() {
                if (!fgets(hdr_buf, MAXLINE, bf)) {
                    io_check(!ferror(bf), "read spill header");
                    return false;
                }
                io_check(fgets(seq_buf, MAXLINE, bf) &&
                         fgets(plus_buf, MAXLINE, bf) &&
                         fgets(qual_buf, MAXLINE, bf), "read complete spill record");
                return true;
            };
            auto at = [&](long off) {
                io_check(fseek(bf, off, SEEK_SET) == 0, "seek spill record");
                io_check(read_record(), "read indexed spill record");
            };
            std::vector<Cand> cands;
            while (true) {
                long off = ftell(bf);
                io_check(off >= 0, "tell spill position");
                if (!read_record()) break;
                int sl = strlen(seq_buf) - 1;
                uint64_t h = 1469598103934665603ULL;
                for (const char* p = hdr_buf + 1; *p && *p != ' ' &&
                     *p != '\t' && *p != '\n'; ++p) {
                    h ^= (unsigned char)*p;
                    h *= 1099511628211ULL;
                }
                cands.push_back({score_read(qual_buf, sl, window_size), h, off, sl});
            }
            std::sort(cands.begin(), cands.end(), [&](const Cand& a, const Cand& c) {
                if (a.score != c.score) return a.score > c.score;
                if (a.idhash != c.idhash) return a.idhash < c.idhash;
                // Rare ID ties/hash collisions: compare full records, not spill
                // offsets. Identical records are interchangeable.
                at(a.off);
                std::string left = std::string(hdr_buf) + seq_buf + plus_buf + qual_buf;
                at(c.off);
                return left < std::string(hdr_buf) + seq_buf + plus_buf + qual_buf;
            });
            for (const Cand& cd : cands) {
                if (accepted_bases == target_bases) break;
                at(cd.off);
                if (dedupe) {
                    const char* end = hdr_buf + 1;
                    while (*end && *end != ' ' && *end != '\t' && *end != '\n') ++end;
                    if (!seen_ids.emplace(hdr_buf + 1, end - hdr_buf - 1).second) {
                        ++dup_reads;
                        continue;
                    }
                }
                if (cd.len > target_bases - accepted_bases) continue;
                write_bytes(hdr_buf, strlen(hdr_buf));
                write_bytes(seq_buf, strlen(seq_buf));
                write_bytes(plus_buf, strlen(plus_buf));
                write_bytes(qual_buf, strlen(qual_buf));
                ++accepted_reads;
                accepted_bases += cd.len;
            }
        }
    }

    if (gz_out) io_check(gzclose(gz_out) == Z_OK, "close gzip output");
    else if (out != stdout) io_check(fclose(out) == 0, "close output");
    else io_check(fflush(out) == 0, "flush stdout");

    free(hdr_buf); free(seq_buf); free(plus_buf); free(qual_buf);

    // Report
    fprintf(stderr, "[fastq_filter] Total reads:  %lld\n", total_reads);
    if (dup_reads)
        fprintf(stderr, "[fastq_filter] Dupes removed: %lld\n", dup_reads);
    if (short_reads)
        fprintf(stderr, "[fastq_filter] Below min_len: %lld\n", short_reads);
    if (qual_reads)
        fprintf(stderr, "[fastq_filter] Below min_q:   %lld\n", qual_reads);
    fprintf(stderr, "[fastq_filter] Accepted:     %lld reads, %lld bases (%.1f Gbp)\n",
            accepted_reads, accepted_bases, accepted_bases / 1e9);
    if (filtering) {
        if (filter_by_bases) {
            fprintf(stderr, "[fastq_filter] Target was:   %lld bases (%.1f Gbp)\n",
                    target_bases, target_bases / 1e9);
            double overshoot = (double)(accepted_bases - target_bases) / target_bases * 100;
            fprintf(stderr, "[fastq_filter] Overshoot:    %+.1f%%\n", overshoot);
        }
        if (filter_by_pct) {
            double actual_pct = n_scored > 0 ? (double)accepted_reads / n_scored * 100 : 0;
            fprintf(stderr, "[fastq_filter] Target pct:   %.1f%% of reads (actual: %.1f%%)\n",
                    keep_percent, actual_pct);
        }
        fprintf(stderr, "[fastq_filter] Final thresh: %.2f (mean=%.2f)\n",
                threshold, score_mean);

        // ASCII histogram: score distribution with kept/tossed stacked bars
        // Collapse 1000 bins into ~40 display bins
        static const int DISP_BINS = 40;
        int bins_per_disp = HIST_BINS / DISP_BINS;
        long long disp_total[DISP_BINS], disp_kept[DISP_BINS];
        long long max_count = 0;
        // Find range with data
        int first_bin = DISP_BINS, last_bin = 0;
        for (int d = 0; d < DISP_BINS; d++) {
            disp_total[d] = 0;
            disp_kept[d] = 0;
            for (int b = d * bins_per_disp; b < (d + 1) * bins_per_disp; b++) {
                disp_total[d] += score_hist[b];
                disp_kept[d] += score_hist_kept[b];
            }
            if (disp_total[d] > 0) {
                if (d < first_bin) first_bin = d;
                if (d > last_bin) last_bin = d;
                if (disp_total[d] > max_count) max_count = disp_total[d];
            }
        }
        if (max_count > 0 && last_bin >= first_bin) {
            static const int BAR_WIDTH = 50;
            fprintf(stderr, "\n[fastq_filter] Score distribution (# kept, . tossed):\n");
            for (int d = first_bin; d <= last_bin; d++) {
                double score_lo = (double)d * bins_per_disp / HIST_BINS * 100.0;
                int bar_total = (int)((double)disp_total[d] / max_count * BAR_WIDTH);
                int bar_kept = (int)((double)disp_kept[d] / max_count * BAR_WIDTH);
                if (disp_total[d] > 0 && bar_total == 0) bar_total = 1;
                if (disp_kept[d] > 0 && bar_kept == 0) bar_kept = 1;
                int bar_tossed = bar_total - bar_kept;
                if (bar_tossed < 0) bar_tossed = 0;

                fprintf(stderr, "  %5.1f |", score_lo);
                for (int i = 0; i < bar_kept; i++) fputc('#', stderr);
                for (int i = 0; i < bar_tossed; i++) fputc('.', stderr);
                fprintf(stderr, " %lld/%lld\n", disp_kept[d], disp_total[d]);
            }
            fprintf(stderr, "         ");
            for (int i = 0; i < BAR_WIDTH + 1; i++) fputc('-', stderr);
            fprintf(stderr, "\n         # = kept, . = tossed\n");
        }
    }

    return 0;
}

int main(int argc, char** argv) {
    // Turn a closed downstream pipe into a checked write failure so RAII can
    // remove spill files before returning a nonzero exit status.
    std::signal(SIGPIPE, SIG_IGN);
    try {
        return run_main(argc, argv);
    } catch (const std::exception& error) {
        fprintf(stderr, "[fastq_filter] ERROR: %s\n", error.what());
        return 1;
    }
}
