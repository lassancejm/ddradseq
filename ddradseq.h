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
#include "khash.h"

/* Set buffer size */
#define BUFLEN 0x20000

/* Define for forward and reverse reads */
#define FORWARD 1
#define REVERSE 2

/* Define operating modes */
enum _runmode_t { PARSE, TRIMEND, PAIR, ERROR };
typedef enum _runmode_t runmode_t;

/* Define data structures */

/* Define command line parameter data structure */
typedef struct cmdparam
{
	bool is_reverse;
	bool default_dir;
	char *parentdir;
	char *outdir;
	char *filename;
	char *csvfile;
	int num_threads;
} CMD;

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
	char *poutfile;
	char *pbuffer;
	size_t pcurr_bytes;
	khash_t(barcode) *b;
} POOL;

/* Second-level hash */
KHASH_MAP_INIT_STR(pool, POOL*)

/* Top-level hash */
KHASH_MAP_INIT_STR(pool_hash, khash_t(pool)*)

/* Function prototypes */

extern CMD* parse_cmdline(int argc, char *argv[], char *mode);

extern int check_directories(CMD *cp, khash_t(pool_hash) *h);

extern khash_t(pool_hash)* read_csv(CMD *cp);

extern void free_db(CMD *cp, khash_t(pool_hash) *h);

extern int parse_fastq(CMD *cp, khash_t(pool_hash) *h);

extern int parse_buffer(CMD *cp, char *buff, const size_t nl, khash_t(pool_hash) *h);

extern char* clean_buffer(char *buff, size_t *nl);

extern size_t reset_buffer(char *buff, const char *r);

extern size_t count_lines(const char *buff);

extern int flush_buffer(BARCODE *bc);

extern int flush_pbuffer(POOL *pl);

extern int print_buffer(char *buff, const size_t nl);

extern int free_cmdline(CMD *cp);

extern int levenshtein(char *s1, char *s2);

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
