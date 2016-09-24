/* file: ddradseq.h
 * description: Header for the ddradseq program
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#ifndef DDRADSEQ_H
#define DDRADEQ_H

#include <stdlib.h>
#include <string.h>
#include "khash.h"

/* Set buffer size */
#define BUFLEN 0x4000

/* Define for forward and reverse reads */
#define FORWARD 1
#define REVERSE 2

/* Define operating modes */
enum mode_t { PARSE, TRIMEND, PAIR, ERROR };

/* Define data structures */

/* Define command line parameter data structure */
typedef struct cmdparam
{
	int default_dir;
	int num_threads;
	char *parentdir;
	char *outdir;
	char *filename;
	char *csvfile;
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
	char *dirpath;
	khash_t(barcode) *b;
} POOL;

/* Second-level hash */
KHASH_MAP_INIT_STR(pool, POOL*)

/* Top-level hash */
KHASH_MAP_INIT_STR(pool_hash, khash_t(pool)*)

/* Function prototypes */

extern CMD* parse_cmdline(int argc, char *argv[], char *mode);

extern int check_directories(CMD *cp, khash_t(pool_hash) *h);

extern khash_t(pool_hash)* read_csv(char *filename, char *outpath);

extern void free_db(khash_t(pool_hash) *h);

extern int parse_fastq(char *filename, khash_t(pool_hash) *h);

extern int parse_buffer(char *buff, const size_t nl, khash_t(pool_hash) *h);

extern char* clean_buffer(char *buff, size_t *nl);

extern size_t reset_buffer(char *buff, const char *r);

extern size_t count_lines(const char *buff);

extern int flush_buffer(BARCODE *bc);

extern int print_buffer(char *buff, const size_t nl);

extern int free_cmdline(CMD *cp);

extern int levenshtein(char *s1, char *s2);

static inline mode_t find_mode(const char *m)
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
