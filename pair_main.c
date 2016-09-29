/* file: pair_main.c
 * description: Entry point for the pair modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ddradseq.h"
#include "khash.h"

void pair_usage(void);

int
pair_main(int argc, char *argv[])
{
	char *pch = NULL;
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
	f = traverse_dirtree(cp->outdir, "parse", &nfiles);

	for (i = 0; i < nfiles; i += 2)
	{
		khash_t(fastq) *h = NULL;
		char *ffor = NULL;
		char *frev = NULL;
		size_t spn = 0;

		/* Construct output file names */
		ffor = malloc(strlen(f[i]) + 1u);
		assert(ffor != NULL);
		frev = malloc(strlen(f[i + 1]) + 1u);
		assert(frev != NULL);
		strcpy(ffor, f[i]);
		strcpy(frev, f[i + 1]);
		pch = strstr(ffor, "parse");
		strncpy(pch, "pairs", 5);
		pch = strstr(frev, "parse");
		strncpy(pch, "pairs", 5);
		spn = strcspn(ffor, ".");
		if (strncmp(ffor, frev, spn) != 0)
		{
			fprintf(stderr, "Error pairing files: %s and %s\n", ffor, frev);
			exit(1);
		}

		/* Read forward fastQ file into hash table */
		h = fastq_to_db(f[i]);

		/* Align mated pairs and write to output file*/
		pair_mates(f[i+1], h, ffor, frev);

		/* Free allocated memory */
		free(ffor);
		free(frev);
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

void
pair_usage(void)
{
	fputs("Usage : ddradseq pair [OPTIONS] [PARENTDIR]\n\n", stderr);
	fputs("Aligns mated pairs in all fastQ files existing in a directory tree.\n\n", stderr);
	fputs("Mandatory arguments to long options are mandatory for short options too.\n", stderr);
	fputs(" -h, --help			 Display this help message\n\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
