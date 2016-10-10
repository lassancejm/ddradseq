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

int pair_main(CMD *cp)
{
	char *pch = NULL;
	char **f = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int n = 0;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Get list of all files */
	f = traverse_dirtree(cp->outdir, "parse", &n);
	if (!f)
		return 1;

	for (i = 0; i < n; i += 2)
	{
		khash_t(fastq) *h = NULL;
		char *ffor = NULL;
		char *frev = NULL;
		size_t spn = 0;
		size_t strl = 0;

		/* Construct output file names */
		strl = strlen(f[i]);
		ffor = malloc(strl + 1u);
		if (UNLIKELY(!ffor))
		{
			logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
		strl = strlen(f[i + 1]);
		frev = malloc(strl + 1u);
		if (UNLIKELY(!frev))
		{
			logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
		strcpy(ffor, f[i]);
		strcpy(frev, f[i + 1]);
		pch = strstr(ffor, "parse");
		if (!pch)
			return 1;
		strncpy(pch, "pairs", 5);
		pch = strstr(frev, "parse");
		if (!pch)
			return 1;
		strncpy(pch, "pairs", 5);

		/* Double check that files are mates */
		spn = strcspn(ffor, ".");
		ret = strncmp(ffor, frev, spn);
		if (ret)
		{
			logerror("%s:%d Pairing files \'%s\' and \'%s\' failed.\n", __func__,
			         __LINE__, ffor, frev);
			return 1;
		}

		/* Read forward fastQ file into hash table */
		h = fastq_to_db(f[i]);
		if (!h)
			return 1;

		/* Print informational update to log file */
		fprintf(lf, "[ddradseq: %s] INFO -- Attempting to pair files \'%s\' and \'%s\'.\n",
		        timestr, ffor, frev);

		/* Align mated pairs and write to output file*/
		ret = pair_mates(f[i + 1], h, ffor, frev);
		if (ret) return 1;

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
	for (i = 0; i < n; i++)
		free(f[i]);
	free(f);

	return 0;
}
