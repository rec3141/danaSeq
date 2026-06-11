/*
 * fastq_filter — Single-pass streaming FASTQ dedup + quality/length filter.
 *
 * Replaces: cat *.fastq.gz | dedup | filtlong -t TARGET
 *
 * Single pass, pipe-friendly. Reads are scored and immediately accepted or
 * rejected using a dynamically adjusted threshold that converges on
 * target_bases (biased to undershoot for memory safety).
 *
 * Scoring follows filtlong's algorithm (without global z-score normalization):
 *   length_score  = 100 * len / (len + 5000)
 *   mean_quality  = mean(1 - 10^(-Q/10)) * 100
 *   window_quality = min sliding-window avg * 100
 *   final = sqrt(length_score * mean_quality) * (0.5 + 0.5 * window_ratio)
 *
 * Input: a custom raw-inflate reader (GzReader) decompresses gzip/BGZF, spans
 * concatenated members (e.g. `cat a.fastq.gz b.fastq.gz`), and RESYNCS past
 * trailing garbage between members rather than silently stopping. (zlib's
 * gzgets stops at the first trailing-garbage byte and can't distinguish that
 * from EOF, which silently truncated cat'd BGZF streams.) Genuine mid-member
 * corruption is reported and yields a non-zero exit.
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
#include <vector>
#include <unordered_set>
#include <getopt.h>
#include <zlib.h>
#include <sys/stat.h>

// Precomputed Phred+33 -> P(correct) lookup
static double PHRED_LUT[256];

static void init_phred_lut() {
    for (int i = 0; i < 256; i++) PHRED_LUT[i] = 0.0;
    for (int q = 0; q < 94; q++) {
        PHRED_LUT[q + 33] = 1.0 - pow(10.0, -q / 10.0);
    }
}

// ============================================================================
// GzReader — robust streaming gzip/BGZF decompressor.
//
// Decompresses a sequence of inputs (files, or "-" for stdin) using raw
// inflate(). Handles the things zlib's gz* API botches on a piped, cat'd
// BGZF stream:
//   - concatenated gzip members (incl. empty BGZF EOF-marker members) are
//     spanned by inflateReset()-ing at each Z_STREAM_END;
//   - trailing garbage between/after members is skipped by scanning forward
//     to the next gzip magic (1f 8b 08) instead of treating it as EOF;
//   - a real inflate error mid-member sets `error` (caller exits non-zero)
//     instead of silently truncating.
// Exposes getline() with the same contract as zlib's gzgets (returns the line
// including its trailing '\n', NUL-terminated; nullptr at end of all input).
// ============================================================================
class GzReader {
    std::vector<const char*> files;
    size_t file_idx = 0;
    FILE* fp = nullptr;
    z_stream strm;
    bool strm_init = false;
    bool all_done = false;

    static const size_t INBUF  = 1u << 20;  // 1 MB compressed read chunk
    static const size_t OUTBUF = 1u << 20;  // 1 MB decompressed window
    unsigned char* inbuf;
    char* outbuf;
    size_t out_pos = 0, out_end = 0;        // consumed / filled in outbuf

    bool open_next_source() {
        if (fp && fp != stdin) fclose(fp);
        fp = nullptr;
        while (file_idx < files.size()) {
            const char* f = files[file_idx++];
            if (strcmp(f, "-") == 0) { fp = stdin; }
            else {
                fp = fopen(f, "rb");
                if (!fp) { fprintf(stderr, "[fastq_filter] Cannot open %s\n", f); continue; }
            }
            if (!strm_init) { memset(&strm, 0, sizeof(strm)); inflateInit2(&strm, 15 + 16); strm_init = true; }
            else            { inflateReset2(&strm, 15 + 16); }
            strm.next_in = inbuf; strm.avail_in = 0;
            return true;
        }
        all_done = true;
        return false;
    }

    // After a Z_STREAM_END, position the input at the next gzip member, skipping
    // any trailing garbage. Returns false if no further member in this source.
    bool start_next_member() {
        while (true) {
            // Need at least 3 bytes to test for the gzip magic.
            if (strm.avail_in < 3) {
                memmove(inbuf, strm.next_in, strm.avail_in);
                size_t got = fread(inbuf + strm.avail_in, 1, INBUF - strm.avail_in, fp);
                strm.next_in = inbuf; strm.avail_in += got;
                if (got == 0) return false;   // source drained
            }
            const unsigned char* p = strm.next_in;
            size_t n = strm.avail_in, i = 0;
            bool found = false;
            for (; i + 3 <= n; i++) {
                if (p[i] == 0x1f && p[i+1] == 0x8b && p[i+2] == 0x08) { found = true; break; }
            }
            if (found) {
                if (i > 0) resyncs++;          // skipped i bytes of trailing garbage
                strm.next_in += i; strm.avail_in -= i;
                inflateReset2(&strm, 15 + 16);
                return true;
            }
            // No magic yet; keep a 2-byte tail (a split magic) and refill.
            size_t keep = (n >= 2) ? 2 : n;
            memmove(inbuf, strm.next_in + n - keep, keep);
            size_t got = fread(inbuf + keep, 1, INBUF - keep, fp);
            strm.next_in = inbuf; strm.avail_in = keep + got;
            if (got == 0) return false;
        }
    }

    // Fill outbuf with as much decompressed data as fits, spanning gzip members
    // (BGZF emits one member per ~64 KB block) so the caller refills rarely and
    // the per-member transition cost is amortised. Returns bytes produced.
    size_t fill() {
        strm.next_out = (unsigned char*)outbuf; strm.avail_out = OUTBUF;
        while (!all_done && strm.avail_out > 0) {
            if (strm.avail_in == 0) {
                size_t got = fread(inbuf, 1, INBUF, fp);
                if (got == 0) { if (!open_next_source()) break; continue; }
                strm.next_in = inbuf; strm.avail_in = got;
            }
            int ret = inflate(&strm, Z_NO_FLUSH);
            if (ret == Z_STREAM_END) {
                if (!start_next_member() && !open_next_source()) all_done = true;
                continue;                      // keep packing the next member in
            }
            if (ret == Z_OK) continue;         // progress; loop refills in/out as needed
            if (ret == Z_BUF_ERROR) {
                if (strm.avail_out == 0) break;        // outbuf full -> return it
                if (strm.avail_in == 0) continue;      // needs more input -> refill at top
                break;                                  // no progress possible
            }
            // Z_DATA_ERROR / Z_NEED_DICT / Z_MEM_ERROR: corruption inside a member.
            fprintf(stderr, "[fastq_filter] decompression error (%d): %s\n",
                    ret, strm.msg ? strm.msg : "?");
            error = true;
            if (!open_next_source()) all_done = true;
        }
        return OUTBUF - strm.avail_out;
    }

public:
    bool error = false;
    long long resyncs = 0;                     // # trailing-garbage skips

    explicit GzReader(const std::vector<const char*>& in) : files(in) {
        inbuf  = (unsigned char*)malloc(INBUF);
        outbuf = (char*)malloc(OUTBUF);
        open_next_source();
    }
    ~GzReader() {
        if (strm_init) inflateEnd(&strm);
        if (fp && fp != stdin) fclose(fp);
        free(inbuf); free(outbuf);
    }

    // gzgets-compatible: copies one line (incl. trailing '\n') into buf,
    // NUL-terminated. Returns buf, or nullptr once all input is consumed.
    char* getline(char* buf, int max) {
        int w = 0;
        while (true) {
            if (out_pos < out_end) {
                void* nl = memchr(outbuf + out_pos, '\n', out_end - out_pos);
                size_t avail = nl ? (size_t)((char*)nl - (outbuf + out_pos) + 1)
                                  : (out_end - out_pos);
                size_t cpy = avail;
                if ((int)cpy > max - 1 - w) cpy = max - 1 - w;   // clamp to buffer
                memcpy(buf + w, outbuf + out_pos, cpy);
                w += cpy; out_pos += avail;                      // consume full line
                if (nl) { buf[w] = '\0'; return buf; }
            } else {
                size_t got = fill();
                if (got == 0) { if (w == 0) return nullptr; buf[w] = '\0'; return buf; }
                out_pos = 0; out_end = got;
            }
        }
    }
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
        "\nReads from stdin if no files given.\n"
        "Deduplication is ON by default (unique to fastq_filter).\n", prog);
}

int main(int argc, char** argv) {
    init_phred_lut();

    long long target_bases = 0;
    double keep_percent = 0.0;  // 0 = disabled
    int min_length = 0;
    double min_mean_q = 0.0;
    double min_window_q = 0.0;
    int window_size = 250;
    bool dedupe = true;
    const char* output_path = nullptr;

    static struct option long_opts[] = {
        {"target_bases",  required_argument, 0, 't'},
        {"keep_percent",  required_argument, 0, 'p'},
        {"min_length",    required_argument, 0, 'm'},
        {"min_mean_q",    required_argument, 0, 'Q'},
        {"min_window_q",  required_argument, 0, 'W'},
        {"window_size",   required_argument, 0, 'w'},
        {"no_dedupe",     no_argument,       0, 'D'},
        {"output",        required_argument, 0, 'o'},
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
            case 'h': usage(argv[0]); return 0;
            default:  usage(argv[0]); return 1;
        }
    }

    bool filter_by_bases = target_bases > 0;
    bool filter_by_pct = keep_percent > 0.0 && keep_percent < 100.0;
    bool filter_by_qual = min_mean_q > 0.0 || min_window_q > 0.0;
    bool filtering = filter_by_bases || filter_by_pct;

    double pct_z = 0.0;
    if (filter_by_pct) {
        pct_z = inv_normal(1.0 - keep_percent / 100.0);
    }
    (void)pct_z;

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
    long long score_hist[HIST_BINS];
    long long score_hist_kept[HIST_BINS];
    memset(score_hist, 0, sizeof(score_hist));
    memset(score_hist_kept, 0, sizeof(score_hist_kept));

    // Line buffers
    const int MAXLINE = 1 << 20; // 1 MB max line
    char* hdr_buf  = (char*)malloc(MAXLINE);
    char* seq_buf  = (char*)malloc(MAXLINE);
    char* plus_buf = (char*)malloc(MAXLINE);
    char* qual_buf = (char*)malloc(MAXLINE);

    auto write_record = [&](const char* hdr, int hlen,
                            const char* seq, int slen,
                            const char* qual, int qlen) {
        if (gz_out) {
            gzwrite(gz_out, hdr, hlen);
            gzwrite(gz_out, "\n", 1);
            gzwrite(gz_out, seq, slen);
            gzwrite(gz_out, "\n+\n", 3);
            gzwrite(gz_out, qual, qlen);
            gzwrite(gz_out, "\n", 1);
        } else {
            fwrite(hdr, 1, hlen, out);
            fputc('\n', out);
            fwrite(seq, 1, slen, out);
            fputs("\n+\n", out);
            fwrite(qual, 1, qlen, out);
            fputc('\n', out);
        }
    };

    GzReader reader(inputs);

    while (true) {
        if (!reader.getline(hdr_buf,  MAXLINE)) break;
        if (!reader.getline(seq_buf,  MAXLINE)) break;
        if (!reader.getline(plus_buf, MAXLINE)) break;
        if (!reader.getline(qual_buf, MAXLINE)) break;

        total_reads++;

        // Strip newlines
        int hlen = strlen(hdr_buf);
        if (hlen > 0 && hdr_buf[hlen-1] == '\n') hdr_buf[--hlen] = '\0';
        int slen = strlen(seq_buf);
        if (slen > 0 && seq_buf[slen-1] == '\n') seq_buf[--slen] = '\0';
        int qlen = strlen(qual_buf);
        if (qlen > 0 && qual_buf[qlen-1] == '\n') qual_buf[--qlen] = '\0';

        // Dedup by read ID
        if (dedupe) {
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

        // Update threshold
        if (n_scored == 10 || (n_scored > 10 && n_scored % 100 == 0)) {
            double score_std = (n_scored > 1) ? sqrt(score_m2 / n_scored) : 0;
            if (score_std > 0) {
                if (filter_by_pct) {
                    double target_frac = keep_percent / 100.0;
                    double actual_frac = (n_scored > 0)
                        ? (double)accepted_reads / n_scored : 1.0;
                    double error = actual_frac - target_frac;
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
                    double remaining_need = (double)target_bases - accepted_bases;
                    double remaining_input = (double)scored_bases;
                    if (remaining_input < 1) remaining_input = 1;

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

    if (gz_out) gzclose(gz_out);
    else if (out != stdout) fclose(out);

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
    if (reader.resyncs > 0)
        fprintf(stderr, "[fastq_filter] WARNING: skipped %lld trailing-garbage region(s) "
                        "between gzip members (corrupt input file(s))\n", reader.resyncs);
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

        // ASCII histogram
        static const int DISP_BINS = 40;
        int bins_per_disp = HIST_BINS / DISP_BINS;
        long long disp_total[DISP_BINS], disp_kept[DISP_BINS];
        long long max_count = 0;
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

    // Genuine mid-member decompression error => non-zero exit (never silently
    // truncate + succeed). Trailing-garbage resyncs are warned about, not fatal.
    return reader.error ? 1 : 0;
}
