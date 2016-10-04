/* file: parse_main.c
 * description: Entry point for the parse modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include "khash.h"
#include "ddradseq.h"

/* Function prototypes */
void parse_deallocate_mem(khash_t(pool_hash)*, khash_t(mates)*, char**, unsigned int, char*, char*);

int
parse_main(CMD *cp)
{
	char **f = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int nfiles = 0;
	khash_t(pool_hash) *h = NULL;
	khash_t(mates) *m = NULL;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Read CSV database into memory */
	h = read_csv(cp);
	if (h == NULL)
	{
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to read CSV database "
		        "into memory.\n", timestr, __func__, __LINE__);
		return 1;
	}

	/* Check for write permissions on parent of output directory */
	if ((ret = check_directories(cp, h)) != 0)
	{
		parse_deallocate_mem(h, NULL, NULL, 0, NULL, NULL);
		return 1;
	}

	/* Initialize hash for mate pair information */
	if ((m = kh_init(mates)) == NULL)
	{
		parse_deallocate_mem(h, NULL, NULL, 0, NULL, NULL);
		return 1;
	}

	/* Get list of all files */
	if ((f = traverse_dirtree(cp->parentdir, NULL, &nfiles)) == NULL)
	{
		parse_deallocate_mem(h, m, NULL, 0, NULL, NULL);
		return 1;
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
			parse_deallocate_mem(h, m, f, nfiles, NULL, NULL);
			return 1;
		}
		if ((frev = malloc(strlen(f[i + 1]) + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			parse_deallocate_mem(h, m, f, nfiles, ffor, NULL);
			return 1;
		}
		strcpy(ffor, f[i]);
		strcpy(frev, f[i + 1]);
		spn = strcspn(ffor, ".");
		if (strncmp(ffor, frev, spn) != 0)
		{
			fprintf(stderr, "ERROR: \'%s\' and \'%s\' do not appear to be mate pairs.\n", ffor, frev);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Files \'%s\' and \'%s\' do not appear to be mate pairs.\n",
			        timestr, __func__, __LINE__, ffor, frev);
			parse_deallocate_mem(h, m, f, nfiles, ffor, frev);
			return 1;
		}

		/* Print informational update to log file */
		fprintf(lf, "[ddradseq: %s] INFO -- Attempting to align sequences in "
		        "\'%s\' and \'%s\'.\n", timestr, ffor, frev);

		/* Read the forward fastQ input file */
		if ((ret = parse_fastq(FORWARD, ffor, h, m, cp->dist)) != 0)
		{
			parse_deallocate_mem(h, m, f, nfiles, ffor, frev);
			return 1;
		}
	
		/* Read the reverse fastQ input file */
		if ((ret = parse_fastq(REVERSE, frev, h, m, cp->dist)) != 0)
		{
			parse_deallocate_mem(h, m, f, nfiles, ffor, frev);
			return 1;
		}
		parse_deallocate_mem(NULL, NULL, NULL, nfiles, ffor, frev);
	}

	/* Deallocate memory */
	parse_deallocate_mem(h, m, f, nfiles, NULL, NULL);

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log */
	fprintf(lf, "[ddradseq: %s] INFO -- Parse step of pipeline is complete.\n",
	        timestr);

	return 0;
}

void
parse_deallocate_mem(khash_t(pool_hash) *h, khash_t(mates) *m, char **f, unsigned int n, char *r, char *s)
{
	if (f)
	{
		unsigned int i = 0;
		for (i = 0; i < n; i++) free(f[i]);
		free(f);
	}
	if (h) free_db(h);
	if (m) kh_destroy(mates, m);
	if (r) free(r);
	if (s) free(s);
}
