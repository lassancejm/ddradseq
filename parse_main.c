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
void deallocate_mem(khash_t(pool_hash)*, khash_t(mates)*);

int
parse_main(CMD *cp)
{
	int ret = 0;
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
		deallocate_mem(h, NULL);
		return 1;
	}

	/* Initialize hash for mate pair information */
	if ((m = kh_init(mates)) == NULL)
	{
		deallocate_mem(cp, h, NULL);
		return 1;
	}

	/* Read the forward fastQ input file */
	if ((ret = parse_fastq(FORWARD, cp->forfastq, h, m, cp->dist)) != 0)
	{
		deallocate_mem(h, m);
		return 1;
	}

	/* Read the reverse fastQ input file */
	if ((ret = parse_fastq(REVERSE, cp->revfastq, h, m, cp->dist)) != 0)
	{
		deallocate_mem(h, m);
		return 1;
	}

	/* Deallocate memory */
	deallocate_mem(h, m);

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log */
	fprintf(lf, "[ddradseq: %s] INFO -- Parse step of pipeline is complete.\n",
	        timestr);

	return 0;
}

void
deallocate_mem(CMD *cp, khash_t(pool_hash) *h, khash_t(mates) *m)
{
	if (h) free_db(h);
	if (m) kh_destroy(mates, m);
}
