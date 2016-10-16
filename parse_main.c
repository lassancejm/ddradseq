/* file: parse_main.c
 * description: Entry point for the parse modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

int parse_main(CMD *cp)
{
	char **f = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int n = 0;
	khash_t(pool_hash) *h = NULL;
	khash_t(mates) *m = NULL;

	/* Read CSV database into memory */
	h = read_csv(cp);
	if (!h)
	{
		logerror("%s:%d Failed to read CSV database into memory.\n",
			     __func__, __LINE__);
		return 1;
	}

	/* Check for write permissions on parent of output directory */
	ret = check_directories(cp, h);
	if (ret)
		return 1;

	/* Initialize hash for mate pair information */
	m = kh_init(mates);
	if (!m)
		return 1;

	/* Get list of all files */
	f = traverse_dirtree(cp->parentdir, NULL, &n);
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
		spn = strcspn(ffor, ".");
		ret = strncmp(ffor, frev, spn);
		if (ret)
		{
			logerror("%s:%d Files \'%s\' and \'%s\' do not appear to be mate-"
			         "pairs.\n", __func__, __LINE__, ffor, frev);
			return 1;
		}

		/* Update time string */
		get_timestr(&timestr[0]);

		/* Print informational update to log file */
		fprintf(lf, "[ddradseq: %s] INFO -- Attempting to align sequences in "
		        "\'%s\' and \'%s\'.\n", timestr, ffor, frev);

		/* Read the forward fastQ input file */
		ret = parse_fastq(FORWARD, ffor, h, m, cp->dist);
		if (ret)
			return 1;

		/* Read the reverse fastQ input file */
		ret = parse_fastq(REVERSE, frev, h, m, cp->dist);
		if (ret)
			return 1;
		free(ffor);
		free(frev);
	}

	/* Deallocate memory from the heap */
	for (i = 0; i < n; i++)
		free(f[i]);
	free(f);
	free_db(h);
	free_matedb(m);

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log */
	fprintf(lf, "[ddradseq: %s] INFO -- Parse step of pipeline is complete.\n",
	        timestr);

	return 0;
}
