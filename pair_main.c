/* file: pair_main.c
 * description: Entry point for the pair modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"
#include "khash.h"

void pair_main_deallocate(char**, char*, char*, unsigned int);

int
pair_main(CMD *cp)
{
	char *pch = NULL;
	char **f = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int nfiles = 0;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Get list of all files */
	if ((f = traverse_dirtree(cp->outdir, "parse", &nfiles)) == NULL)
		return 1;

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
			pair_main_deallocate(f, NULL, NULL, nfiles);
			return 1;
		}
		strl = strlen(f[i + 1]);
		if ((frev = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			pair_main_deallocate(f, ffor, NULL, nfiles);
			return 1;
		}
		strcpy(ffor, f[i]);
		strcpy(frev, f[i + 1]);
		pch = strstr(ffor, "parse");
		strncpy(pch, "pairs", 5);
		pch = strstr(frev, "parse");
		strncpy(pch, "pairs", 5);

		/* Double check that files are mates */
		spn = strcspn(ffor, ".");
		if ((ret = strncmp(ffor, frev, spn)) != 0)
		{
			fprintf(stderr, "ERROR: Pairing files \'%s\' and \'%s\' failed.\n", ffor, frev);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Pairing files \'%s\' and \'%s\' failed.\n",
			        timestr, __func__, __LINE__, ffor, frev);
			pair_main_deallocate(f, ffor, frev, nfiles);
			return 1;
		}

		/* Read forward fastQ file into hash table */
		if ((h = fastq_to_db(f[i])) == NULL)
		{
			pair_main_deallocate(f, ffor, frev, nfiles);
			return 1;
		}

		/* Print informational update to log file */
		fprintf(lf, "[ddradseq: %s] INFO -- Attempting to pair files \'%s\' and \'%s\'.\n",
		        timestr, ffor, frev);

		/* Align mated pairs and write to output file*/
		if ((ret = pair_mates(f[i + 1], h, ffor, frev)) != 0)
		{
			free_pairdb(h);
			pair_main_deallocate(f, ffor, frev, nfiles);
			return 1;
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
	pair_main_deallocate(f, NULL, NULL, nfiles);

	return 0;
}


void
pair_main_deallocate(char **f, char *ff, char *rf, unsigned int n)
{
	if (f)
	{
		unsigned int i = 0;
		for (i = 0; i < n; i++) free(f[i]);
		free(f);
	}
	if (ff) free(ff);
	if (rf) free(rf);
}
