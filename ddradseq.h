/* file: ddradseq.h
 * description: Header for the ddradseq program
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#ifndef DDRADSEQ_H
#define DDRADEQ_H

#if defined __GNUC__ && !defined __GNUC_STDC_INLINE__ && !defined __GNUC_GNU_INLINE__
#define __GNUC_GNU_INLINE__ 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
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

/* File I/O buffer size */
#define BUFLEN 0x20000

/* Maximum line length to read from input file */
#define MAX_LINE_LENGTH 400

/* Number of lines in individual parse buffers */
#define BSIZE 4000

/* Length of terminal output directory name */
#define DNAME_LENGTH 5

/* Define for forward and reverse reads */
#define FORWARD 1
#define REVERSE 2

/* Length of data format YYYY-DD-MM */
#define DATELEN 20

/* Define alignment status bit positions */
#define KSW_XBYTE  0x10000
#define KSW_XSTOP  0x20000
#define KSW_XSUBO  0x40000
#define KSW_XSTART 0x80000


/******************************************************
 * Data structure defintions
 ******************************************************/

/* Data structure to hold command line parameters */
typedef struct cmdparam
{
	bool across;
	char *parent_indir;
	char *parent_outdir;
	char *outdir;
	char *csvfile;
	char *mode;
	char *glob;
	int dist;
	int score;
	int gapo;
	int gape;
	int nthreads;
	FILE *lf;
} CMD;

/* Data structure to hold a single fastQ entry */
typedef struct _fastq_
{
	char *id;
	char *seq;
	char *qual;
} FASTQ;

/* Data structure to hold the results of a local sequence alignment */
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

/* Data structure to hold align query parameters */
typedef struct _kswq_t
{
    int qlen;
    int slen;
    unsigned char shift;
    unsigned char mdiff;
    unsigned char max;
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


/******************************************************
 * Function prototypes
 ******************************************************/

/******************************************************
 * Parsing functions
 ******************************************************/

/* int parse_fastq(const CMD *cp, const int orient, const char *filename, khash_t(pool_hash) *h,
 *                 khash_t(mates) *m)
 * Parses a fastQ file by index sequence
 * Arguments:
 * cp -- Pointer to command line data structure (read-only)
 * orient -- Orientation of reads in fastQ file (read-only)
 * filename -- Pointer to string holding fastQ input file name (read-only)
 * h -- Pointer to pool_hash hash table with parsing database
 * m -- Pointer to mate information hash table
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int parse_fastq(const CMD *cp, const int orient, const char *filename,
                       khash_t(pool_hash) *h, khash_t(mates) *m);


/* int parse_forwardbuffer(const CMD *cp, char *buff, const size_t nl, khash_t(pool_hash) *h,
 *                         khash_t(mates) *m)
 * Parses forward fastQ entries in the buffer
 * Arguments:
 * cp -- Pointer to command line data structure (read-only)
 * buff -- Pointer to string holding the buffer
 * nl -- Number of lines in the buffer (read-only)
 * h -- Pointer to pool_hash hash table with parsing database (read-only)
 * m -- Pointer to mate information hash table
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int parse_forwardbuffer(const CMD *cp, char *buff, const size_t nl, const khash_t(pool_hash) *h,
                               khash_t(mates) *m);


/* int parse_reversebuffer(const CMD *cp, char *buff, const size_t nl, khash_t(pool_hash) *h,
                           khash_t(mates) *m)
 * Parses reverse fastQ entries in the buffer
 * Arguments:
 * cp -- Pointer to command line data structure (read-only)
 * buff -- Pointer to string holding the buffer
 * nl -- Number of lines in the buffer (read-only)
 * h -- Pointer to pool_hash hash table with parsing database (read-only)
 * m -- Pointer to mate information hash table (read-only)
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int parse_reversebuffer(const CMD *cp, char *buff, const size_t nl, const khash_t(pool_hash) *h,
                               const khash_t(mates) *m);


/******************************************************
 * Sequence pairing functions
 ******************************************************/

/* int pair_mates(const char *filename, const khash_t(fastq) *h, const char *ffor, const char *frev, FILE *lf)
 * Pairs mates in two fastQ files
 * Arguments:
 * filename -- Pointer to string for input forward fastQ (read-only)
 * h -- Pointer to hash table to hold forward sequences (read only)
 * ffor -- Pointer to string with forward output file name (read only)
 * frev -- Pointer to string with reverse output file name (read only)
 * lf -- Pointer to log file stream
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int pair_mates(const char *filename, const khash_t(fastq) *h, const char *ffor, const char *frev, FILE *lf);


/******************************************************
 * Trimend functions
 ******************************************************/

/* int align_mates(const CMD *cp, const char *fin, const char *rin, const char *fout, const char *rout)
 * Align mates in two fastQ files and trim 3' end of reverse sequences
 * Arguments:
 * cp -- Pointer to command line data structure (read-only)
 * fin -- Pointer to string with forward input file name (read-only)
 * rin -- Pointer to string with reverse input file name (read-only)
 * fout -- Pointer to string with forward output file name (read-only)
 * rout -- Pointer to string with reverse output file name (read-only)
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int align_mates(const CMD *cp, const char *fin, const char *rin, const char *fout, const char *rout);


/******************************************************
 * UI functions
 ******************************************************/

/* CMD *get_cmdline(int argc, char *argv[])
 * Reads command line parameters into CMD data structure
 * Arguments:
 * argc -- Number of command line arguments
 * argv -- Pointer to array of strings holding command line arguments
 * Returns:
 * Pointer to command line data structure on success or NULL on failure
 */

extern CMD *get_cmdline(int argc, char *argv[]);


/******************************************************
 * I/O management functions
 ******************************************************/

/* khash_t(pool_hash) *read_csv(const CMD *cp)
 * Reads CSV database file into parsing hash database
 * Arguments:
 * cp -- Pointer to command line data structure (read-only)
 * Returns:
 * Pointer to pool_hash hash table on success or NULL on failure
 */

extern khash_t(pool_hash) *read_csv(const CMD *cp);


/* int check_csv(const CMD *cp)
 * Check the integrity of the input CSV file
 * Arguments:
 * cp -- Pointer to the command line parameter data structure (read-only)
 * Returns:
 * Zero on pass, non-zero on fail integrity test
 */

extern int check_csv(const CMD *cp);


/* khash_t(fastq)* fastq_to_db(const char *filename, FILE *lf)
 * Populates a fastQ database from fastQ input file
 * Arguments:
 * filename -- Pointer to string holding input fastQ file name (read-only)
 * lf -- Pointer to log file stream
 * Returns:
 * Pointer to fastQ hash table on success or NULL on failure
 */

extern khash_t(fastq) *fastq_to_db(const char *filename, FILE *lf);


/******************************************************
 * File system functions
 ******************************************************/

/* int create_dirtree(const CMD *cp, const khash_t(pool_hash) *h)
 * Creates and checks output directory tree
 * Arguments:
 * cp -- Pointer to command line data structure (read-only)
 * h -- Pointer to pool_hash hash table with parsing database (read-only)
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int create_dirtree(const CMD *cp, const khash_t(pool_hash) *h);


/* unsigned int traverse_dirtree(const CMD *cp, const char *caller, char ***flist)
 * Produces a sorted list of all fastQ files in the input directory tree
 * Arguments:
 * cp -- Pointer to command line data structure (read-only)
 * caller -- Pointer to string identifying calling function (read-only)
 * flist -- Pointer to array of input file names
 * Returns:
 * The number of input files found in the input directory tree
 */

extern unsigned int traverse_dirtree(const CMD *cp, const char *caller, char ***flist);


/******************************************************
 * Buffer management functions
 ******************************************************/

/* char *clean_buffer(char *buff, size_t *nl)
 * Limits the buffer to hold only entire fastQ entries
 * Arguments:
 * buff -- Pointer to the string holding the buffer
 * nl -- Number of new line characters in the buffer
 * Returns:
 * Pointer to new end of buffer
 */

extern char *clean_buffer(char *buff, size_t *nl);


/* size_t reset_buffer(char *buff, const char *r)
 * Resets the buffer for next read block
 * Arguments:
 * buff -- Pointer to the string holding the buffer
 * r -- Pointer to new end of buffer (read-only)
 * Returns:
 * The number of leftover characters in the old buffer
 */

extern size_t reset_buffer(char *buff, const char *r);


/* size_t reset_buffer(char *buff, const char *r)
 * Counts newline characters in buffer
 * Arguments:
 * buff -- Pointer to the string holding the buffer (read-only)
 * Returns:
 * The number of newline characters in the buffer
 */

extern size_t count_lines(const char *buff);


/* int flush_buffer(int orient, BARCODE *bc)
 * Dumps a full buffer to file
 * Arguments:
 * orient -- Orientation of reads in the buffer
 * bc -- Pointer to BARCODE data structure
 * lf -- Pointer to log file stream
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int flush_buffer(int orient, BARCODE *bc, FILE *lf);


/******************************************************
 * Memory management functions
 ******************************************************/

/* int destroy_cmdline(CMD *cp)
 * Destroy command line parameter data structure
 * Arguments:
 * cp -- Pointer to command line data structure
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int destroy_cmdline(CMD *cp);


/* int free_db(khash_t(pool_hash) *h)
 * Deallocates memory used by CSV database
 * Arguments:
 * h -- Pointer to pool_hash hash table with parsing database
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int free_db(khash_t(pool_hash) *h);


/* int free_pairdb (khash_t(fastq) *h)
 * Deallocates memory used by forward fastQ database
 * Arguments:
 * h -- Pointer to hash table to hold forward sequences
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int free_pairdb(khash_t(fastq) *h);


/* int free_matedb(khash_t(mates) *m)
 * Deallocates memory used by mate pair database
 * Arguments:
 * m -- Pointer to mate information hash table
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int free_matedb(khash_t(mates) *m);


/******************************************************
 * Alignment functions
 ******************************************************/

/* ALIGN_RESULT local_align(int qlen, char *query, int tlen, char *target, const char *mat,
                            int gapo, int gape, int xtra, FILE *lf)
 * Calculates the local sequence alignment by Smith-Waterman algorithm
 * Arguments:
 * qlen -- Length of query sequence
 * query -- Pointer string holding query sequence
 * tlen -- Length of the target sequence
 * target -- Pointer to string holding target sequence
 * mat -- Scoring matrix in a one-dimension array (read-only)
 * gapo -- Gap penalty
 * gape -- Gap extension penalty
 * xtra -- Status variable
 * lf -- Pointer to log file stream
 * Returns:
 * ALIGN_RESULT data structure
 */

extern ALIGN_RESULT local_align(int qlen, char *query, int tlen, char *target, const char *mat,
                                int gapo, int gape, int xtra, FILE *lf);


/* char* revcom(const char *s, FILE *lf)
 * Reverse complement a DNA string with full IUPAC alphabet
 * Arguments:
 * s -- Pointer to string to be reverse-complemented (read-only)
 * lf -- Pointer to log file stream
 * Returns:
 * Pointer to newly allocated reverse-complemented string on success or NULL on failure
 */

extern char *revcom(const char *s, FILE *lf);


/******************************************************
 * Edit distance functions
 ******************************************************/

/* int levenshtein(const char *s1, const char *s2)
 * Calculates the Levenshtein distance between two strings
 * Arguments:
 * s1 -- Pointer to the first string (read-only)
 * s2 -- Pointer to the second string (read-only)
 * Returns:
 * Edit distance between the two strings
 */

extern int levenshtein(const char *s1, const char *s2);


/******************************************************
 * Log file functions
 ******************************************************/

/* int log_init(CMD *cp)
 * Initialize the ddradseq log file
 * Arguments:
 * cp -- Pointer to command line data structure (read-only)
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int log_init(CMD *cp);


/* void loginfo(FILE *lf, const char *format)
 * Write informational message to log file
 */

extern void loginfo(FILE *lf, const char *format, ...);


/* void logwarn(FILE *lf, const char *format)
 * Write warning message to log file
 */

extern void logwarn(FILE *lf, const char *format, ...);


/* void logerror(FILE *lf, const char *format)
 * Report error to both logfile and stderr
 */

extern void logerror(FILE *lf, const char *format, ...);


/* int get_timestr(char *s);
 * Fill in current time for log file reporting
 * Arguments:
 * s -- Pointer to time string
 * Returns:
 * Zero on success and non-zero on failure
 */

extern int get_timestr(char *s);

/******************************************************
 * Error reporting functions
 ******************************************************/

/* void error(const char *format)
 * Report error to stderr
 */

extern void error(const char *format, ...);


/******************************************************
 * Inline utility functions
 ******************************************************/

static inline int string_equal(const char *a, const char *b)
{
	return (strcmp(a, b) == 0);
}

#endif
