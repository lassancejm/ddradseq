/* file: pair_main.c
 * description: Entry point for the pair modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include "ddradseq.h"
#include "khash.h"

void pair_usage(void);
static int compare(const void * a, const void * b);

int
pair_main(int argc, char *argv[])
{
	char **f = NULL;
	unsigned int i = 0;
	unsigned int nfiles = 0;
	CMD *cp = NULL;

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "pair")) == NULL)
	{
		pair_usage();
		return 1;
	}

	/* Get list of all files */	
	f = traverse_dirtree(cp->outdir, &nfiles);

	/* Sort file list */
	qsort (f, nfiles, sizeof(const char *), compare);

	for (i = 0; i < nfiles; i += 2)
	{
		khash_t(fastq) *h = NULL;

		/* Read forward fastQ file into hash table */
		h = fastq_to_db(f[i]);
	
		/* Align mated pairs and write to output file*/
		pair_mates(f[i+1], h);
		free_pairdb(h);
	}

	/* Deallocate memory */
	for (i = 0; i < nfiles; i++)
		free(f[i]);
	free(f);
	if (cp)
		free_cmdline(cp);

	return 0;
}

static int
compare(const void * a, const void * b)
{
    /* The pointers point to offsets into "array", so we need to
       dereference them to get at the strings. */

    return strcmp(*(const char **) a, *(const char **) b);
}

void
pair_usage(void)
{
	fputs("Usage : ddradseq pair [OPTIONS] [FASTQ.R1] [FASTQ.R2]\n\n", stderr);
	fputs("Aligns mated pairs in two fastQ files.\n\n", stderr);
	fputs("Mandatory arguments to long options are mandatory for short options too.\n", stderr);
	fputs(" -o, --out=DIR        Parent directory to write output files    [default: same as input fastQ]\n", stderr);
	fputs(" -n, --threads=INT    Number of threads for concurrency         [default: 1]\n", stderr);
	fputs(" -h, --help           Display this help message\n\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
