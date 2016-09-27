/* file: free_db.c
 * description: Function to free memory used by CSV database
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include "khash.h"
#include "ddradseq.h"

void
free_db (khash_t(pool_hash) *h)
{
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khash_t(pool) *p = NULL;
	khash_t(barcode) *b = NULL;
	POOL *pl = NULL;
	BARCODE *bc = NULL;

	if (h == NULL) return;

	for (i = kh_begin(h); i != kh_end(h); i++)
	{
		if (kh_exist(h, i))
		{
			p = kh_value(h, i);
			for (j = kh_begin(p); j != kh_end(p); j++)
			{
				if (kh_exist(p, j))
				{
					pl = kh_value(p, j);
					free(pl->poolID);
					free(pl->poolpath);
					b = pl->b;
					for (k = kh_begin(b); k != kh_end(b); k++)
					{
						if (kh_exist(b, k))
						{
							bc = kh_value(b, k);
							free(bc->smplID);
							free(bc->outfile);
							free(bc->buffer);
						}
					}
					kh_destroy(barcode, b);
				}
			}
			kh_destroy(pool, p);
		}
	}
	kh_destroy(pool_hash, h);
}
