/* file: ddradseq.h
 * description: Header for the ddradseq program
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#ifndef DDRADSEQ_H
#define DDRADEQ_H

#if defined __GNUC__ && !defined __GNUC_STDC_INLINE__ && !defined __GNUC_GNU_INLINE__
#define __GNUC_GNU_INLINE__ 1
#endif

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <emmintrin.h>
#include "khash.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/* Set buffer size */
#define BUFLEN 0x20000

/* Define for forward and reverse reads */
#define FORWARD 1
#define REVERSE 2

#define DATELEN 11

/* Define operating modes */
enum _runmode_t { PARSE, TRIMEND, PAIR, ERROR };
typedef enum _runmode_t runmode_t;

/* Define data structures */

/* Define command line parameter data structure */
typedef struct cmdparam
{
	bool default_dir;
	bool trim_barcode;
	char *parentdir;
	char *outdir;
	char *forfastq;
	char *revfastq;
	char *csvfile;
	int num_threads;
} CMD;

typedef struct _fastq_
{
	char *id;
	char *seq;
	char *qual;
} FASTQ;

typedef struct _ksqr_t
{
	int score;
	int target_begin;
	int target_end;
	int query_begin;
	int query_end;
	int score2;
	int target_end2;
} ALIGN_RESULT;

typedef struct _kswq_t
{
    int qlen;
    int slen;
    unsigned char shift;
    unsigned char mdiff;
    unsigned char max;
    unsigned char size;
    __m128i *qp;
    __m128i *H0;
    __m128i *H1;
    __m128i *E;
    __m128i *Hmax;
} ALIGN_QUERY;

/* Barcode-level data structure */
typedef struct _barcode_
{
	char *smplID;
	char *outfile;
	char *buffer;
	size_t curr_bytes;
} BARCODE;

/* Third-level hash */
KHASH_MAP_INIT_STR(barcode, BARCODE*)

/* Pool-level data structure */
typedef struct _pool_
{
	char *poolID;
	char *poolpath;
	size_t barcode_length;
	khash_t(barcode) *b;
} POOL;

/* Second-level hash */
KHASH_MAP_INIT_STR(pool, POOL*)

/* Top-level hash */
KHASH_MAP_INIT_STR(pool_hash, khash_t(pool)*)

/* Hash to hold fastQ entries */
KHASH_MAP_INIT_STR(fastq, FASTQ*)

/* Hash to hold mate information */
KHASH_MAP_INIT_STR(mates, char*)

/* Function prototypes */

extern CMD* parse_cmdline(int argc, char *argv[], char *mode);

extern int check_directories(CMD *cp, khash_t(pool_hash) *h);

extern khash_t(pool_hash)* read_csv(CMD *cp);

extern void free_db(khash_t(pool_hash) *h);

extern void free_pairdb (khash_t(fastq) *h);

extern int parse_fastq(int orient, char *filename, khash_t(pool_hash) *h,
                       khash_t(mates) *m, bool trim_barcode);

extern int parse_forwardbuffer(char *buff, const size_t nl, khash_t(pool_hash) *h,
                               khash_t(mates) *m, bool trim_barcode);

extern int parse_reversebuffer(char *buff, const size_t nl, khash_t(pool_hash) *h,
                               khash_t(mates) *m);

extern char* clean_buffer(char *buff, size_t *nl);

extern size_t reset_buffer(char *buff, const char *r);

extern size_t count_lines(const char *buff);

extern char** traverse_dirtree(const char *dirpath, int *x);

extern int flush_buffer(int orient, BARCODE *bc);

extern int free_cmdline(CMD *cp);

extern int levenshtein(char *s1, char *s2);

extern khash_t(fastq)* fastq_to_db(char *filename);

extern char* revcom(char *input_string);

extern int pair_mates(char *filename, khash_t(fastq) *h, char *ffor, char *frev);

extern ALIGN_RESULT local_align(int qlen, unsigned char *query, int tlen,
                                unsigned char *target, int m, const char *mat,
                                int gapo, int gape, int xtra, ALIGN_QUERY **qry);

static inline runmode_t find_mode(const char *m)
{
	int ret;
	if ((ret = strcmp(m, "parse")) == 0)
		return PARSE;
	else if ((ret = strcmp(m, "trimend")) == 0)
		return TRIMEND;
	else if ((ret = strcmp(m, "pair")) == 0)
		return PAIR;
	else
		return ERROR;
}

#endif
