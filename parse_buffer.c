/* file: parse_buffer.c
 * description: Function to parse fastQ entries in the buffer
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

int
parse_buffer(char *buff, const size_t nl, khash_t(pool_hash) *h)
{
	char *q = buff;
	char *r = NULL;
	char *copy = NULL;
	char *idline = NULL;
    char seps[] = ": ";
    char *tok = NULL;
    char *flowcell_ID = NULL;
    char *index_sequence = NULL;
    char *barcode_sequence = NULL;
    char *dna_sequence = NULL;
    char *qual_sequence = NULL;
    int read = 0;
	int ret = 0;
    int z = 0;
    int orient = 0;
	size_t add_bytes = 0;
    size_t strl = 0;
	size_t l = 0;
	size_t ll = 0;
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khint_t m = 0;
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

				/* Get lane number, tile number, x, and y coordinate */
				for (z = 0; z < 4; z++)
				{
					tok = strtok_r(NULL, seps, &r);
					assert(tok != NULL);
				}

				/* Get read number */
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				read = atoi(tok);
				/* Determine orientation of reads in file */
				/* Ideally this approach is wrong and the */
				/* program should be able to handle fastQs */
				/* with mixed read orientations */
				if (read == FORWARD)
				{
					if (orient == 0) orient = FORWARD;
					else if (orient == REVERSE) orient = -1;
				}
				if (read == REVERSE)
				{
					if (orient == 0) orient = REVERSE;
					else if (orient == FORWARD) orient = -1;
				}

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
				if (read == FORWARD)
				{
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
						for (m = kh_begin(b); m != kh_end(b); m++)
						{
							if (kh_exist(b, m))
							{
								int d = 0;
								if ((d = levenshtein((char*)kh_key(b, m), barcode_sequence)) < 2)
								{
									bc = kh_value(b, m);
									break;
								}
							}
						}
					}
				}
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
				if ((bc->curr_bytes + add_bytes) >= BUFLEN)
				{
					if ((ret = flush_buffer(bc)) != 0)
					{
						fputs("Problem writing buffer to file.\n", stderr);
						abort();
					}
				}
				bc->curr_bytes += add_bytes;
				char *t = malloc(add_bytes + 1u);
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
	return orient;
}
