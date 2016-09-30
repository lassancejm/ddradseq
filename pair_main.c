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
#include "ddradseq.h"
#include "khash.h"

void pair_usage(void);
void pair_main_deallocate(CMD*, char**, char*, char*, unsigned int);

int
pair_main(int argc, char *argv[])
{
	char *pch = NULL;
	char **f = NULL;
	int ret = EXIT_SUCCESS;
	unsigned int i = 0;
	unsigned int nfiles = 0;
	CMD *cp = NULL;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "pair")) == NULL)
	{
		pair_usage();
		return EXIT_FAILURE;
	}

	/* Get list of all files */
	if ((f = traverse_dirtree(cp->outdir, "parse", &nfiles)) == NULL)
	{
		free_cmdline(cp);
		return EXIT_FAILURE;
	}

	for (i = 0; i < nfiles; i += 2)
	{
		khash_t(fastq) *h = NULL;
		char *ffor = NULL;
		char *frev = NULL;
		size_t spn = 0;
		size_t strl = 0;

		/* Construct output file names */
		strl = strlen(f[i]);
		if ((ffor = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			pair_main_deallocate(cp, f, NULL, NULL, nfiles);
			return EXIT_FAILURE;
		}
		strl = strlen(f[i + 1]);
		if ((frev = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			pair_main_deallocate(cp, f, ffor, NULL, nfiles);
			return EXIT_FAILURE;
		}
		strcpy(ffor, f[i]);
		strcpy(frev, f[i + 1]);
		pch = strstr(ffor, "parse");
		strncpy(pch, "pairs", 5);
		pch = strstr(frev, "parse");
		strncpy(pch, "pairs", 5);

		/* Double check that files are mates */
		spn = strcspn(ffor, ".");
		if (strncmp(ffor, frev, spn) != 0)
		{
			fprintf(stderr, "ERROR: Pairing files \'%s\' and \'%s\' failed.\n", ffor, frev);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Pairing files \'%s\' and \'%s\' failed.\n",
			        timestr, __func__, __LINE__, ffor, frev);
			pair_main_deallocate(cp, f, ffor, frev, nfiles);
			return EXIT_FAILURE;
		}

		/* Read forward fastQ file into hash table */
		if ((h = fastq_to_db(f[i])) == NULL)
		{
			pair_main_deallocate(cp, f, ffor, frev, nfiles);
			return EXIT_FAILURE;			
		}

		/* Print informational update to log file */
		fprintf(lf, "[ddradseq: %s] INFO -- Attempting to pair files \'%s\' and \'%s\'.\n",
		        timestr, ffor, frev);

		/* Align mated pairs and write to output file*/
		if ((ret = pair_mates(f[i + 1], h, ffor, frev)) == EXIT_FAILURE)
		{
			free_pairdb(h);
			pair_main_deallocate(cp, f, ffor, frev, nfiles);
			return EXIT_FAILURE;			
		}

		/* Free allocated memory */
		free(ffor);
		free(frev);
		free_pairdb(h);
	}

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Print informational message to logfile */
	fprintf(lf, "[ddradseq: %s] INFO -- Done pairing all fastQ files in \'%s\'.\n",
	        timestr, cp->outdir);

	/* Deallocate memory */
	pair_main_deallocate(cp, f, NULL, NULL, nfiles);

	return EXIT_SUCCESS;
}


void
pair_main_deallocate(CMD *cp, char **f, char *ff, char *rf, unsigned int n)
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
pair_usage(void)
{
	fputs("Usage : ddradseq pair [OPTIONS] [PARENTDIR]\n\n", stderr);
	fputs("Aligns mated pairs in all fastQ files existing in a directory tree.\n\n", stderr);
	fputs("Mandatory arguments to long options are mandatory for short options too.\n", stderr);
	fputs(" -h, --help			 Display this help message\n\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
