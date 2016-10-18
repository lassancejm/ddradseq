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

int pair_main(const CMD *cp)
{
	char *pch = NULL;
	char **f = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int n = 0;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Get list of all files */
	if (string_equal(cp->mode, "pair"))
		f = traverse_dirtree(cp->parent_indir, "parse", &n);
	else
		f = traverse_dirtree(cp->outdir, "parse", &n);
	if (!f)
		return 1;

	for (i = 0; i < n; i += 2)
	{
		khash_t(fastq) *h = NULL;
		char *ffor = NULL;
		char *frev = NULL;
		size_t spn = 0;

		/* Construct output file names */
		ffor = strdup(f[i]);
		if (UNLIKELY(!ffor))
		{
			logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
		frev = strdup(f[i+1]);
		if (UNLIKELY(!frev))
		{
			logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
		pch = strstr(ffor, "parse");
		if (!pch)
			return 1;
		strncpy(pch, "pairs", DNAME_LENGTH);
		pch = strstr(frev, "parse");
		if (!pch)
			return 1;
		strncpy(pch, "pairs", DNAME_LENGTH);

		/* Double-check that files are mates */
		spn = strcspn(ffor, ".");
		ret = strncmp(ffor, frev, spn);
		if (ret)
		{
			logerror("%s:%d Files \'%s\' and \'%s\' do not appear to be mate-"
				     "pairs.\n", __func__, __LINE__, ffor, frev);
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
	if (string_equal(cp->mode, "pair"))
		fprintf(lf, "[ddradseq: %s] INFO -- Done pairing all fastQ files in \'%s\'.\n",
		        timestr, cp->parent_indir);
	else
		fprintf(lf, "[ddradseq: %s] INFO -- Done pairing all fastQ files in \'%s\'.\n",
		        timestr, cp->outdir);

	/* Deallocate memory */
	for (i = 0; i < n; i++)
		free(f[i]);
	free(f);

	return 0;
}
