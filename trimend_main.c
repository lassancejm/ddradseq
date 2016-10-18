/* file: trimend_main.c
 * description: Entry point for the trimend function
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

int trimend_main(const CMD *cp)
{
	char *pch = NULL;
	char **f = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int n = 0;

	/* Get time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log file */
	fprintf(lf, "[ddradseq: %s] INFO -- Beginning to trim 3\' end of reverse "
	        "sequences in \'%s\'.\n", timestr, cp->outdir);

	/* Get list of all files */
	if (string_equal(cp->mode, "trimend"))
		f = traverse_dirtree(cp->parent_indir, "pairs", &n);
	else
		f = traverse_dirtree(cp->outdir, "pairs", &n);
	if (!f)
		return 1;

	for (i = 0; i < n; i += 2)
	{
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
		pch = strstr(ffor, "pairs");
		if (!pch)
			return 1;
		strncpy(pch, "final", DNAME_LENGTH);
		pch = strstr(frev, "pairs");
		if (!pch)
			return 1;
		strncpy(pch, "final", DNAME_LENGTH);

		/* Double-check that files are mates */
		spn = strcspn(ffor, ".");
		ret = strncmp(ffor, frev, spn);
		if (ret)
		{
			logerror("%s:%d Files \'%s\' and \'%s\' do not appear to be mate-"
				     "pairs.\n", __func__, __LINE__, ffor, frev);
			return 1;
		}

		/* Print informational update to log file */
		fprintf(lf, "[ddradseq: %s] INFO -- Attempting to align sequences in "
		        "\'%s\' and \'%s\'.\n", timestr, ffor, frev);

		/* Align mated pairs and write to output file*/
		ret = align_mates(cp, f[i], f[i+1], ffor, frev);
		if (ret)
			return 1;

		/* Free allocated memory */
		free(ffor);
		free(frev);
	}

	/* Get time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log file */
	if (string_equal(cp->mode, "trimend"))
	fprintf(lf, "[ddradseq: %s] INFO -- Done trimming 3\' end of reverse "
		        "sequences in \'%s\'.\n", timestr, cp->parent_indir);
	else
		fprintf(lf, "[ddradseq: %s] INFO -- Done trimming 3\' end of reverse "
		        "sequences in \'%s\'.\n", timestr, cp->outdir);

	/* Deallocate memory */
	for (i = 0; i < n; i++)
		free(f[i]);
	free(f);

	return 0;
}
