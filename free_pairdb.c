/* file: free_pairdb.c
 * description: Function to free memory used by pairs database
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include "khash.h"
#include "ddradseq.h"

int
free_pairdb (khash_t(fastq) *h)
{
	khint_t i = 0;
	FASTQ *e = NULL;

	if (h == NULL) return 1;

	for (i = kh_begin(h); i != kh_end(h); i++)
	{
		if (kh_exist(h, i))
		{
			e = kh_value(h, i);
			free(e->id);
			free(e->seq);
			free(e->qual);
		}
	}
	kh_destroy(fastq, h);
	return 0;
}
