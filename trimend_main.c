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
void trimend_main_deallocate(char**, char*, char*, unsigned int);

int
trimend_main(CMD *cp)
{
	char *pch = NULL;
	char **f = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int nfiles = 0;

	/* Get time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log file */
	fprintf(lf, "[ddradseq: %s] INFO -- Beginning to trim 3\' end of reverse "
	        "sequences in \'%s\'.\n", timestr, cp->outdir);

	/* Get list of all files */
	if ((f = traverse_dirtree(cp->outdir, "pairs", &nfiles)) == NULL)
		return 1;

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
			trimend_main_deallocate(f, NULL, NULL, nfiles);
			return 1;
		}
		if ((frev = malloc(strlen(f[i + 1]) + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			trimend_main_deallocate(f, ffor, NULL, nfiles);
			return 1;
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
			trimend_main_deallocate(f, ffor, frev, nfiles);
			return 1;
		}

		/* Print informational update to log file */
		fprintf(lf, "[ddradseq: %s] INFO -- Attempting to align sequences in "
		        "\'%s\' and \'%s\'.\n", timestr, ffor, frev);

		/* Align mated pairs and write to output file*/
		if ((ret = align_mates(f[i], f[i + 1], ffor, frev)) != 0)
		{
			trimend_main_deallocate(f, ffor, frev, nfiles);
			return 1;			
		}

		/* Free allocated memory */
		free(ffor);
		free(frev);
	}

	/* Get time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log file */
	fprintf(lf, "[ddradseq: %s] INFO -- Done trimming 3\' end of reverse "
	        "sequences in \'%s\'.\n", timestr, cp->outdir);

	/* Deallocate memory */
	trimend_main_deallocate(f, NULL, NULL, nfiles);

	return 0;
}

void
trimend_main_deallocate(char **f, char *ff, char *rf, unsigned int n)
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
