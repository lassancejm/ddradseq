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
#include "khash.h"
#include "ddradseq.h"

#define KEYLEN 31

int
parse_forwardbuffer(char *buff, const size_t nl, khash_t(pool_hash) *h, khash_t(mates) *m, int dist)
{
	char *q = buff;
	char *s = NULL;
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
	unsigned char *skip = NULL;
	int a = 0;
	int ret = 0;
	int z = 0;
	int tile = 0;
	int xpos = 0;
	int ypos = 0;
	size_t add_bytes = 0;
	size_t strl = 0;
	size_t l = 0;
	size_t ll = 0;
	size_t sl = 0;
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khint_t kk = 0;
	khint_t mk = 0;
	khash_t(barcode) *b = NULL;
	khash_t(pool) *p = NULL;
	BARCODE *bc = NULL;
	POOL *pl = NULL;

	if ((skip = malloc(nl)) == NULL)
	{
		fputs("ERROR: cannot allocate memory for skip array.\n", stderr);
		exit(EXIT_FAILURE);
	}
	for (l = 0; l < nl; l++) skip[l] = 0;

	for (l = 0; l < nl; l++)
	{
		ll = strlen(q);
		if (skip[l] == 0)
		{
			switch (l % 4)
			{
				case 0:
					/* Illumina identifier line */
					/* Make a copy of the Illumina identifier line */
					if ((copy = malloc(ll + 1u)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for copy.\n", stderr);
						exit(EXIT_FAILURE);
					}
					if ((idline = malloc(ll + 1u)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for idline.\n", stderr);
						exit(EXIT_FAILURE);
					}
					strcpy(idline, q);
					strcpy(copy, q);

					/* Get instrument name */
					if ((tok = strtok_r(copy, seps, &r)) == NULL)
					{
						fputs("ERROR: strtok_r failed.\n", stderr);
						exit(EXIT_FAILURE);
					}

					/* Get run number */
					if ((tok = strtok_r(NULL, seps, &r)) == NULL)
					{
						fputs("ERROR: strtok_r failed.\n", stderr);
						exit(EXIT_FAILURE);
					}

					/* Get flow cell ID */
					if ((tok = strtok_r(NULL, seps, &r)) == NULL)
					{
						fputs("ERROR: strtok_r failed.\n", stderr);
						exit(EXIT_FAILURE);
					}
					strl = strlen(tok);
					if ((flowcell_ID = malloc(strl + 1u)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for flowcell_ID.\n", stderr);
						exit(EXIT_FAILURE);
					}
					strcpy(flowcell_ID, tok);
					i = kh_get(pool_hash, h, flowcell_ID);
					p = kh_value(h, i);

					/* Get lane number */
					if ((tok = strtok_r(NULL, seps, &r)) == NULL)
					{
						fputs("ERROR: strtok_r failed.\n", stderr);
						exit(EXIT_FAILURE);
					}

					/* Get tile number, x, and y coordinate */
					if ((tok = strtok_r(NULL, seps, &r)) == NULL)
					{
						fputs("ERROR: strtok_r failed.\n", stderr);
						exit(EXIT_FAILURE);
					}
					tile = atoi(tok);
					if ((tok = strtok_r(NULL, seps, &r)) == NULL)
					{
						fputs("ERROR: strtok_r failed.\n", stderr);
						exit(EXIT_FAILURE);
					}
					xpos = atoi(tok);
					if ((tok = strtok_r(NULL, seps, &r)) == NULL)
					{
						fputs("ERROR: strtok_r failed.\n", stderr);
						exit(EXIT_FAILURE);
					}
					ypos = atoi(tok);

					/* Get read orientation, filtered flag and control number */
					for (z = 0; z < 3; z++)
					{
						if ((tok = strtok_r(NULL, seps, &r)) == NULL)
						{
							fputs("ERROR: strtok_r failed.\n", stderr);
							exit(EXIT_FAILURE);
						}
					}

					/* Get the index sequence */
					if ((tok = strtok_r(NULL, seps, &r)) == NULL)
					{
						fputs("ERROR: strtok_r failed.\n", stderr);
						exit(EXIT_FAILURE);
					}
					strl = strlen(tok);
					if ((index_sequence = malloc(strl + 1u)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for index_sequence.\n", stderr);
						exit(EXIT_FAILURE);
					}
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
					/* Grab the barcode before trimming */
					if ((barcode_sequence = malloc(pl->barcode_length + 1u)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for barcode_sequence.\n", stderr);
						exit(EXIT_FAILURE);
					}
					strncpy(barcode_sequence, q, pl->barcode_length);
					barcode_sequence[pl->barcode_length] = '\0';

					/* Trim the barcode */
					sl = ll - pl->barcode_length;
					if ((dna_sequence = malloc(sl + 1u)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for dna_sequence.\n", stderr);
						exit(EXIT_FAILURE);
					}
					s = q;
					s += pl->barcode_length;
					strcpy(dna_sequence, s);

					/* Find the barcode in the database */
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
								if ((d = levenshtein((char*)kh_key(b, kk),
													 barcode_sequence)) <= dist)
								{
									bc = kh_value(b, kk);
									break;
								}
							}
						}
					}
					/* If barcode still not found-- skip sequence */
					if (bc == NULL)
					{
						skip[l + 1] = 1;
						skip[l + 2] = 1;
						free(idline);
						free(dna_sequence);
						free(barcode_sequence);
						break;
					}

					/* Constrcut key for mate pair hash */
					if ((mkey = malloc(KEYLEN)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for mkey.\n", stderr);
						exit(EXIT_FAILURE);
					}
					sprintf(mkey, "%010d%010d%010d", tile, xpos, ypos);
					mk = kh_put(mates, m, mkey, &a);
					if (!a) free(mkey);
					kh_value(m, mk) = strdup(barcode_sequence);
					break;
				case 2:
					/* Quality identifier line */
					break;
				case 3:
					/* Quality sequence line */
					sl = ll - pl->barcode_length;
					if ((qual_sequence = malloc(sl + 1u)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for qual_sequence.\n", stderr);
						exit(EXIT_FAILURE);
					}
					s = q;
					s += pl->barcode_length;
					strcpy(qual_sequence, s);
					add_bytes = strlen(idline) + strlen(dna_sequence) +
								strlen(qual_sequence) + 5u;
					char *t = NULL;
					if ((bc->curr_bytes + add_bytes) >= BUFLEN)
					{
						if ((ret = flush_buffer(FORWARD, bc)) != 0)
						{
							fputs("Problem writing buffer to file.\n", stderr);
							exit(EXIT_FAILURE);
						}
					}
					bc->curr_bytes += add_bytes;
					if ((t = malloc(add_bytes + 1u)) == NULL)
					{
						fputs("ERROR: cannot allocate memory for t.\n", stderr);
						exit(EXIT_FAILURE);
					}
					sprintf (t, "%s\n%s\n+\n%s\n", idline, dna_sequence,
							 qual_sequence);
					strcat(bc->buffer, t);

					/* Free alloc'd memory for fastQ entry */
					free(t);
					free(idline);
					free(dna_sequence);
					free(qual_sequence);
					free(barcode_sequence);
					break;
			}
		}
		q += ll + 1u;
	}
	free(skip);
	return 0;
}
