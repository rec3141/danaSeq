/*
 * tetramer_freqs — Tetranucleotide frequency calculator for metagenomic contigs.
 *
 * Drop-in replacement for tetramer_freqs.py / tetramer_freqs_esom.pl.
 * Outputs headerless TSV: contig_id<TAB>136 RC-collapsed tetramer frequencies.
 *
 * Compile: gcc -O3 -o tetramer_freqs tetramer_freqs.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

/* 136 canonical (reverse-complement-collapsed) tetranucleotides.
 * Generated at startup by build_canonical_kmers(). */
static int N_CANONICAL;
static int canonical[256]; /* kmer index -> canonical index (0..135) */
static char canonical_names[136][5];

/* Encode a base to 0-3, or -1 for ambiguous */
static inline int base_encode(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

/* Encode a 4-mer string to an index 0-255 */
static inline int kmer_encode(const char *s) {
    int v = 0;
    for (int i = 0; i < 4; i++) {
        int b = base_encode(s[i]);
        if (b < 0) return -1;
        v = (v << 2) | b;
    }
    return v;
}

/* Reverse complement of a 4-mer index */
static int kmer_revcomp(int idx) {
    int rc = 0;
    for (int i = 0; i < 4; i++) {
        rc = (rc << 2) | (3 - (idx & 3));
        idx >>= 2;
    }
    return rc;
}

static void index_to_str(int idx, char *out) {
    const char bases[] = "ACGT";
    for (int i = 3; i >= 0; i--) {
        out[i] = bases[idx & 3];
        idx >>= 2;
    }
    out[4] = '\0';
}

static void build_canonical_kmers(void) {
    int seen[256];
    memset(seen, -1, sizeof(seen));
    N_CANONICAL = 0;
    for (int i = 0; i < 256; i++) {
        if (seen[i] >= 0) {
            canonical[i] = seen[i];
            continue;
        }
        int rc = kmer_revcomp(i);
        canonical[i] = N_CANONICAL;
        seen[i] = N_CANONICAL;
        if (rc != i) {
            canonical[rc] = N_CANONICAL;
            seen[rc] = N_CANONICAL;
        }
        index_to_str(i, canonical_names[N_CANONICAL]);
        N_CANONICAL++;
    }
}

/* Dynamic buffer for reading sequences */
typedef struct {
    char *data;
    size_t len;
    size_t cap;
} Buffer;

static void buf_init(Buffer *b) {
    b->cap = 1 << 20; /* 1 MB initial */
    b->data = malloc(b->cap);
    b->len = 0;
}

static void buf_clear(Buffer *b) { b->len = 0; }

static void buf_push(Buffer *b, const char *s, size_t n) {
    while (b->len + n >= b->cap) {
        b->cap *= 2;
        b->data = realloc(b->data, b->cap);
    }
    memcpy(b->data + b->len, s, n);
    b->len += n;
}

static void buf_free(Buffer *b) { free(b->data); }

/* Mask short regions between ambiguous bases with N */
static void mask_ambiguous(char *seq, size_t len, int min_region) {
    /* Find positions of ambiguous bases */
    size_t *ambig = malloc((len + 2) * sizeof(size_t));
    int n_ambig = 0;
    ambig[n_ambig++] = 0;
    for (size_t i = 1; i < len; i++) {
        if (base_encode(seq[i]) < 0)
            ambig[n_ambig++] = i;
    }
    ambig[n_ambig++] = len;

    for (int i = 1; i < n_ambig; i++) {
        size_t span = ambig[i] - ambig[i - 1] + 1;
        if ((int)span < min_region) {
            for (size_t j = ambig[i - 1]; j < ambig[i]; j++)
                seq[j] = 'N';
        }
    }
    free(ambig);
}

/* Count tetranucleotide frequencies in a sequence region */
static void count_tetras(const char *seq, size_t len, double *freqs) {
    int counts[136];
    memset(counts, 0, sizeof(counts));
    int total = 0;

    /* Streaming 4-mer encoding: maintain a rolling 2-bit encoded value */
    int val = 0;
    int valid = 0; /* number of consecutive valid bases */
    for (size_t i = 0; i < len; i++) {
        int b = base_encode(seq[i]);
        if (b < 0) {
            valid = 0;
            val = 0;
            continue;
        }
        val = ((val << 2) | b) & 0xFF;
        valid++;
        if (valid >= 4) {
            counts[canonical[val]]++;
            total++;
        }
    }

    if (total == 0) {
        for (int i = 0; i < N_CANONICAL; i++) freqs[i] = 0.0;
        return;
    }
    double inv = 1.0 / total;
    for (int i = 0; i < N_CANONICAL; i++)
        freqs[i] = counts[i] * inv;
}

static void write_row(FILE *out, const char *id, double *freqs) {
    fputs(id, out);
    for (int i = 0; i < N_CANONICAL; i++) {
        fputc('\t', out);
        fprintf(out, "%.10g", freqs[i]);
    }
    fputc('\n', out);
}

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s -f input.fasta [-min 2500] [-max 5000] [-o output.lrn] [--lengths lengths.tsv]\n"
        "\n"
        "  -f, --fasta FILE    Input FASTA file\n"
        "  -min N              Minimum contig length [2500]\n"
        "  -max N              Window size for splitting [5000]\n"
        "  -o, --output FILE   Output file [stdout]\n"
        "  --lengths FILE      Output contig lengths\n"
        "  -h, --help          Show this help\n", prog);
}

int main(int argc, char **argv) {
    const char *fasta_path = NULL;
    const char *output_path = NULL;
    const char *lengths_path = NULL;
    int min_len = 2500;
    int max_len = 5000;

    static struct option long_opts[] = {
        {"fasta",   required_argument, 0, 'f'},
        {"output",  required_argument, 0, 'o'},
        {"lengths", required_argument, 0, 'L'},
        {"help",    no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    /* Custom parsing to handle -min and -max (non-standard long opts) */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-min") == 0 && i + 1 < argc) {
            min_len = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-max") == 0 && i + 1 < argc) {
            max_len = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            fasta_path = argv[++i];
        } else if (strcmp(argv[i], "--fasta") == 0 && i + 1 < argc) {
            fasta_path = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_path = argv[++i];
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_path = argv[++i];
        } else if (strcmp(argv[i], "--lengths") == 0 && i + 1 < argc) {
            lengths_path = argv[++i];
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        }
    }

    if (!fasta_path) {
        fprintf(stderr, "Error: -f/--fasta is required\n");
        usage(argv[0]);
        return 1;
    }

    build_canonical_kmers();

    FILE *fasta = fopen(fasta_path, "r");
    if (!fasta) {
        fprintf(stderr, "Error: cannot open %s\n", fasta_path);
        return 1;
    }

    FILE *out = output_path ? fopen(output_path, "w") : stdout;
    if (!out) {
        fprintf(stderr, "Error: cannot open %s\n", output_path);
        return 1;
    }
    FILE *len_out = lengths_path ? fopen(lengths_path, "w") : NULL;

    Buffer seq;
    buf_init(&seq);
    char header[4096];
    header[0] = '\0';
    char line[1 << 20]; /* 1 MB line buffer */
    double freqs[136];
    long long n_contigs = 0, n_written = 0;

    /* Process a completed contig */
    #define FLUSH_CONTIG() do { \
        if (header[0] && (int)seq.len >= min_len) { \
            seq.data[seq.len] = '\0'; \
            mask_ambiguous(seq.data, seq.len, 50); \
            /* Split into windows */ \
            if ((int)seq.len < 2 * max_len) { \
                count_tetras(seq.data, seq.len, freqs); \
                write_row(out, header, freqs); \
                if (len_out) fprintf(len_out, "%s\t%zu\n", header, seq.len); \
                n_written++; \
            } else { \
                int n_frags = 0; \
                size_t frag_starts[65536]; \
                size_t frag_lens[65536]; \
                for (size_t pos = 0; pos < seq.len && n_frags < 65536; pos += max_len) { \
                    frag_starts[n_frags] = pos; \
                    size_t flen = (pos + max_len <= seq.len) ? max_len : seq.len - pos; \
                    frag_lens[n_frags] = flen; \
                    n_frags++; \
                } \
                /* Merge short tail into previous fragment */ \
                if (n_frags > 1 && (int)frag_lens[n_frags-1] < max_len) { \
                    frag_lens[n_frags-2] += frag_lens[n_frags-1]; \
                    n_frags--; \
                } \
                for (int fi = 0; fi < n_frags; fi++) { \
                    char frag_id[4200]; \
                    if (n_frags == 1) \
                        snprintf(frag_id, sizeof(frag_id), "%s", header); \
                    else \
                        snprintf(frag_id, sizeof(frag_id), "%s_%d", header, fi + 1); \
                    count_tetras(seq.data + frag_starts[fi], frag_lens[fi], freqs); \
                    write_row(out, frag_id, freqs); \
                    if (len_out) fprintf(len_out, "%s\t%zu\n", frag_id, frag_lens[fi]); \
                    n_written++; \
                } \
            } \
        } \
    } while(0)

    while (fgets(line, sizeof(line), fasta)) {
        size_t llen = strlen(line);
        if (llen > 0 && line[llen - 1] == '\n') line[--llen] = '\0';
        if (llen > 0 && line[llen - 1] == '\r') line[--llen] = '\0';

        if (line[0] == '>') {
            FLUSH_CONTIG();
            n_contigs++;
            /* Extract header: skip '>', take first space-delimited field */
            char *space = strchr(line + 1, ' ');
            if (space) *space = '\0';
            char *tab = strchr(line + 1, '\t');
            if (tab) *tab = '\0';
            strncpy(header, line + 1, sizeof(header) - 1);
            header[sizeof(header) - 1] = '\0';
            buf_clear(&seq);
        } else {
            buf_push(&seq, line, llen);
        }
    }
    FLUSH_CONTIG();

    fclose(fasta);
    if (out != stdout) fclose(out);
    if (len_out) fclose(len_out);
    buf_free(&seq);

    fprintf(stderr, "[INFO] tetramer_freqs: %lld contigs read, %lld written\n",
            n_contigs, n_written);
    return 0;
}
