/* file: parse_forwardbuffer.c
 * description: Function to parse forward fastQ entries in the buffer
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

#define BARCODE_LENGTH 5
#define KEYLEN 31

int
parse_forwardbuffer(char *buff, const size_t nl, khash_t(pool_hash) *h, khash_t(mates) *m)
{
	char *q = buff;
	char *r = NULL;
	char *copy = NULL;
	char *idline = NULL;
    char seps[] = ": ";
    char *tok = NULL;
    char *mkey = NULL;
    char *flowcell_ID = NULL;
    char *index_sequence = NULL;
    char *barcode_sequence = NULL;
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
				barcode_sequence = malloc(BARCODE_LENGTH + 1u);
				assert(barcode_sequence != NULL);
				strncpy(barcode_sequence, dna_sequence, BARCODE_LENGTH);
				barcode_sequence[BARCODE_LENGTH] = '\0';
				k = kh_get(barcode, b, barcode_sequence);
				if (k != kh_end(b))
					bc = kh_value(b, k);
				else
				{
					/* Iterate through all barcode hash keys and */
					/* calculate Levenshtein distance */
					for (kk = kh_begin(b); kk != kh_end(b); kk++)
					{
						if (kh_exist(b, kk))
						{
							int d = 0;
							if ((d = levenshtein((char*)kh_key(b, kk), barcode_sequence)) < 2)
							{
								bc = kh_value(b, kk);
								break;
							}
						}
					}
				}
				/* If barcode still not found */
				if (bc == NULL) abort();
				/* Need to find a way to skip sequences */

				/* Constrcut key for mate pair hash */
				mkey = malloc(KEYLEN);
				assert(mkey != NULL);
				sprintf(mkey, "%010d%010d%010d", tile, xpos, ypos);
				mk = kh_put(mates, m, mkey, &a);
				if (!a)
					free(mkey);
				kh_value(m, mk) = strdup(barcode_sequence);
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
					if ((ret = flush_buffer(FORWARD, bc)) != 0)
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
