/* file: parse_reversebuffer.c
 * description: Function to parse reverse fastQ entries in the buffer
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "khash.h"
#include "ddradseq.h"

#define KEYLEN 31

int
parse_reversebuffer(char *buff, const size_t nl, khash_t(pool_hash) *h, khash_t(mates) *m)
{
	char *q = buff;
	char *r = NULL;
	char *copy = NULL;
	char *idline = NULL;
    char seps[] = ": ";
    char *tok = NULL;
    char *mkey = NULL;
    char *flowcell_ID = NULL;
    char *barcode_sequence = NULL;
    char *index_sequence = NULL;
    char *dna_sequence = NULL;
    char *qual_sequence = NULL;
    int a = 0;
    int read = 0;
	int ret = 0;
    int z = 0;
    int tile = 0;
    int xpos = 0;
    int ypos = 0;
	size_t add_bytes = 0;
    size_t strl = 0;
	size_t l = 0;
	size_t ll = 0;
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khint_t kk = 0;
	khint_t mk = 0;
	khash_t(barcode) *b = NULL;
	khash_t(pool) *p = NULL;
	BARCODE *bc = NULL;
	POOL *pl = NULL;

	for (l = 0; l < nl; l++)
	{
		ll = strlen(q);
		switch (l % 4)
		{
			case 0:
				/* Illumina identifier line */
				/* Make a copy of the Illumina identifier line */
				copy = malloc(ll + 1u);
				assert(copy != NULL);
				idline = malloc(ll + 1u);
				assert(idline != NULL);
				strcpy(idline, q);
				strcpy(copy, q);

				/* Get instrument name */
				tok = strtok_r(copy, seps, &r);
				assert(tok != NULL);

				/* Get run number */
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);

				/* Get flow cell ID */
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				strl = strlen(tok);
				flowcell_ID = malloc(strl + 1u);
				assert(flowcell_ID != NULL);
				strcpy(flowcell_ID, tok);
				i = kh_get(pool_hash, h, flowcell_ID);
				if (i == kh_end(h))
				{
					fputs("Hash lookup failure.\n", stderr);
					abort();
				}
				p = kh_value(h, i);

				/* Get lane number */
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);

				/* Get tile number, x, and y coordinate */
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				tile = atoi(tok);
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				xpos = atoi(tok);
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				ypos = atoi(tok);

				/* Retrieve barcode sequence of mate */
				mkey = malloc(KEYLEN);
				assert(mkey != NULL);
				sprintf(mkey, "%010d%010d%010d", tile, xpos, ypos);
				mk = kh_get(mates, m, mkey);
				if (mk == kh_end(m))
				{
					fputs("Hash lookup failure.\n", stderr);
					abort();					
				}
				barcode_sequence = kh_value(m, mk);
				free(mkey);

				/* Get read number */
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				read = atoi(tok);

				/* Get filtered flag and control number */
				for (z = 0; z < 2; z++)
				{
					tok = strtok_r(NULL, seps, &r);
					assert(tok != NULL);
				}

				/* Get the index sequence */
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				strl = strlen(tok);
				index_sequence = malloc(strl + 1u);
				assert(index_sequence != NULL);
				strcpy(index_sequence, tok);
				j = kh_get(pool, p, index_sequence);
				pl = kh_value(p, j);
				b = pl->b;

				/* Get the barcode entry of read's mate */
				k = kh_get(barcode, b, barcode_sequence);
				if (k != kh_end(b))
					bc = kh_value(b, k);

				/* Free memory */
				free(flowcell_ID);
				free(index_sequence);
				free(copy);
				break;
			case 1:
				/* Sequence line */
				dna_sequence = malloc(ll + 1u);
				assert(dna_sequence != NULL);
				strcpy(dna_sequence, q);
				break;
			case 2:
				/* Quality identifier line */
				break;
			case 3:
				/* Quality sequence line */
				qual_sequence = malloc(ll + 1u);
				assert(qual_sequence != NULL);
				strcpy(qual_sequence, q);

				add_bytes = strlen(idline) + strlen(dna_sequence) +
				            strlen(qual_sequence) + 5u;
				char *t = NULL;
				if ((bc->curr_bytes + add_bytes) >= BUFLEN)
				{
					if ((ret = flush_buffer(REVERSE, bc)) != 0)
					{
						fputs("Problem writing buffer to file.\n", stderr);
						abort();
					}
				}
				bc->curr_bytes += add_bytes;
				t = malloc(add_bytes + 1u);
				sprintf (t, "%s\n%s\n+\n%s\n", idline, dna_sequence, qual_sequence);
				strcat(bc->buffer, t);

				/* Free alloc'd memory for fastQ entry */
				free(t);
				free(idline);
				free(dna_sequence);
				free(qual_sequence);
				free(barcode_sequence);
				break;
		}
		q += ll + 1u;
	}
	return 0;
}