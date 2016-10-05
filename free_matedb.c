/* file: free_matedb.c
 * description: Function to free memory used by mates database
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include "khash.h"
#include "ddradseq.h"

int
free_matedb(khash_t(mates) *m)
{
	const char *key;
	char *v;

	if (m == NULL) return 1;
	kh_foreach(m, key, v, free(v); free((void*)key););
	kh_destroy(mates, m);
	return 0;
}
