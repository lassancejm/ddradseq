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
void trimend_main_deallocate(CMD*, char**, char*, char*, unsigned int);

int
trimend_main(int argc, char *argv[])
{
	char *pch = NULL;
	char **f = NULL;
	int ret = EXIT_SUCCESS;
	unsigned int i = 0;
	unsigned int nfiles = 0;
	CMD *cp = NULL;

	/* Get time string */
	get_timestr(&timestr[0]);

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "trimend")) == NULL)
	{
		trimend_usage();
		return EXIT_FAILURE;
	}

	/* Print informational message to log file */
	fprintf(lf, "[ddradseq: %s] INFO -- Beginning to trim 3\' end of reverse sequences in \'%s\'.\n",
	        timestr, cp->outdir);

	/* Get list of all files */
	if ((f = traverse_dirtree(cp->outdir, "pairs", &nfiles)) == NULL)
	{
		free_cmdline(cp);
		return EXIT_FAILURE;
	}

	for (i = 0; i < nfiles; i += 2)
	{
		char *ffor = NULL;
		char *frev = NULL;
		size_t spn = 0;

		/* Construct output file names */
		if ((ffor = malloc(strlen(f[i]) + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			trimend_main_deallocate(cp, f, NULL, NULL, nfiles);
			return EXIT_FAILURE;
		}
		if ((frev = malloc(strlen(f[i + 1]) + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			trimend_main_deallocate(cp, f, ffor, NULL, nfiles);
			return EXIT_FAILURE;
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
			fprintf(stderr, "ERROR: \'%s\' and \'%s\' do not appear to be mate pairs.\n", ffor, frev);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Files \'%s\' and \'%s\' do not appear to be mate pairs.\n",
			        timestr, __func__, __LINE__, ffor, frev);
			trimend_main_deallocate(cp, f, ffor, frev, nfiles);
			return EXIT_FAILURE;
		}

		/* Print informational update to log file */
		fprintf(lf, "[ddradseq: %s] INFO -- Attempting to align sequences in \'%s\' and \'%s\'.\n",
		        timestr, ffor, frev);

		/* Align mated pairs and write to output file*/
		if ((ret = align_mates(f[i], f[i + 1], ffor, frev)) == EXIT_FAILURE)
		{
			trimend_main_deallocate(cp, f, ffor, frev, nfiles);
			return EXIT_FAILURE;			
		}

		/* Free allocated memory */
		free(ffor);
		free(frev);
	}

	/* Get time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log file */
	fprintf(lf, "[ddradseq: %s] INFO -- Done trimming 3\' end of reverse sequences in \'%s\'.\n",
	        timestr, cp->outdir);

	/* Deallocate memory */
	trimend_main_deallocate(cp, f, NULL, NULL, nfiles);

	return EXIT_SUCCESS;
}

void
trimend_main_deallocate(CMD *cp, char **f, char *ff, char *rf, unsigned int n)
{
	if (f)
	{
		unsigned int i = 0;
		for (i = 0; i < n; i++) free(f[i]);
		free(f);
	}
	if (ff) free(ff);
	if (rf) free(rf);
	if (cp) free_cmdline(cp);	
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
