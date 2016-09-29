/* file: trimend_main.c
 * description: Entry point for the trimend function
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"
#include "khash.h"

/* Function prototypes */
void trimend_usage (void);

int
trimend_main(int argc, char *argv[])
{
	char *pch = NULL;
	char **f = NULL;
	unsigned int i = 0;
	unsigned int nfiles = 0;
	CMD *cp = NULL;

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "trimend")) == NULL)
	{
		trimend_usage();
		exit(EXIT_FAILURE);
	}

	/* Get list of all files */
	f = traverse_dirtree(cp->outdir, "pairs", &nfiles);

	for (i = 0; i < nfiles; i += 2)
	{
		char *ffor = NULL;
		char *frev = NULL;
		size_t spn = 0;

		/* Construct output file names */
		if ((ffor = malloc(strlen(f[i]) + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for ffor.\n", stderr);
			exit(EXIT_FAILURE);
		}
		if ((frev = malloc(strlen(f[i + 1]) + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for frev.\n", stderr);
			exit(EXIT_FAILURE);
		}
		strcpy(ffor, f[i]);
		strcpy(frev, f[i + 1]);
		pch = strstr(ffor, "pairs");
		strncpy(pch, "trime", 5);
		pch = strstr(frev, "pairs");
		strncpy(pch, "trime", 5);
		spn = strcspn(ffor, ".");
		if (strncmp(ffor, frev, spn) != 0)
		{
			fprintf(stderr, "Error pairing files: %s and %s\n", ffor, frev);
			exit(1);
		}

		/* Align mated pairs and write to output file*/
		align_mates(f[i], f[i + 1], ffor, frev);

		/* Free allocated memory */
		free(ffor);
		free(frev);
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
trimend_usage(void)
{
	fputs("Usage : ddradseq trimend [OPTIONS] [PARENTDIR]\n\n", stderr);
	fputs("Trims the 5\' end of reverse reads in all mate pair fastQ files\n", stderr);
	fputs("in the specified directory tree. Mated pairs are aligned and any\n", stderr);
	fputs(" overhand is trimmed.\n\n", stderr);
	fputs("Mandatory arguments to long options are mandatory for short options too.\n", stderr);
	fputs(" -h, --help			 Display this help message\n\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
