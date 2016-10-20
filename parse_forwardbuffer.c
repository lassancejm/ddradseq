/* file: parse_forwardbuffer.c
 * description: Parses forward fastQ entries in the buffer
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "khash.h"
#include "ddradseq.h"

int parse_forwardbuffer(char *buff, const size_t nl, const khash_t(pool_hash) *h,
                        khash_t(mates) *m, const int dist)
{
	char *q = buff;
	char *s = NULL;
	char *r = NULL;
	char *copy = NULL;
	char *idline = NULL;
	char seps[] = ": ";
	char *tok = NULL;
	char *mkey = NULL;
	char *flowcell = NULL;
	char *index_sequence = NULL;
	char *barcode_sequence = NULL;
	char *dna_sequence = NULL;
	char *qual_sequence = NULL;
	unsigned char *skip = NULL;
	int a = 0;
	int x = 0;
	int ret = 0;
	int tile = 0;
	int xpos = 0;
	int ypos = 0;
	size_t add_bytes = 0;
	size_t l = 0;
	size_t ll = 0;
	size_t sl = 0;
	size_t strl = 0;
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khint_t kk = 0;
	khint_t mk = 0;
	khash_t(barcode) *b = NULL;
	khash_t(pool) *p = NULL;
	BARCODE *bc = NULL;
	POOL *pl = NULL;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Indicator variable whether to skip processing a line */
	skip = calloc(1, nl * sizeof(unsigned char));
	if (!skip)
	{
		logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 1;
	}

	/* Iterate through lines in the buffer */
	for (l = 0; l < nl; l++)
	{
		ll = strlen(q);
		if (!skip[l])
		{
			switch (l % 4)
			{
				case 0:
					/* Make a copy of the Illumina identifier line */
					copy = strdup(q);
					if (!copy)
					{
						error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					idline = strdup(q);
					if (!idline)
					{
						error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}

					for (tok = strtok_r(copy, seps, &r), x = 0; tok != NULL; tok = strtok_r(NULL, seps, &r), x++)
					{
						switch (x)
						{
							case 2:
								flowcell = strdup(tok);
								break;
							case 4:
								tile = atoi(tok);
								break;
							case 5:
								xpos = atoi(tok);
								break;
							case 6:
								ypos = atoi(tok);
								break;
							case 10:
								index_sequence = strdup(tok);
								break;
						}
					}

					/* Check parsing */
					if (!flowcell)
					{
						error("%s:%d Illumina ID parsing failure.\n", __func__, __LINE__);
						return 1;
					}
					if (!index_sequence)
					{
						error("%s:%d Illumina ID parsing failure.\n", __func__, __LINE__);
						return 1;
					}

					/* Lookup flow cell identifier */
					i = kh_get(pool_hash, h, flowcell);
					p = kh_value(h, i);

					/* Lookup pool identifier */
					j = kh_get(pool, p, index_sequence);
					pl = kh_value(p, j);
					b = pl->b;

					/* Free memory */
					free(index_sequence);
					free(copy);
					break;
				case 1:
					/* Sequence line */
					/* Grab the barcode before trimming */
					barcode_sequence = malloc(pl->barcode_length + 1u);
					if (!barcode_sequence)
					{
						error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					strncpy(barcode_sequence, q, pl->barcode_length);
					barcode_sequence[pl->barcode_length] = '\0';

					/* Trim the barcode */
					sl = ll - pl->barcode_length;
					dna_sequence = malloc(sl + 1u);
					if (!dna_sequence)
					{
						error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
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
					if (!bc)
					{
						skip[l+1] = 1;
						skip[l+2] = 1;
						free(idline);
						free(dna_sequence);
						free(barcode_sequence);
						break;
					}

					/* Constrcut key for mate pair hash */
					mkey = malloc(KEYLEN);
					if (!mkey)
					{
						error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					strl = strlen(flowcell);
					sprintf(mkey, "%.*s%s%05d%06d%08d", strl >= 11 ? 0 : (int)(11-strl), "000000000000", flowcell, tile, xpos, ypos);
					mk = kh_put(mates, m, mkey, &a);
					if (a)
						kh_value(m, mk) = strdup(barcode_sequence);
					else
						free(mkey);
					free(flowcell);
					break;
				case 2:
					/* Quality identifier line */
					break;
				case 3:
					/* Quality sequence line */
					sl = ll - pl->barcode_length;
					qual_sequence = malloc(sl + 1u);
					if (!qual_sequence)
					{
						error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					s = q;
					s += pl->barcode_length;
					strcpy(qual_sequence, s);
					add_bytes = strlen(idline) + strlen(dna_sequence) +
								strlen(qual_sequence) + 5u;
					char *t = NULL;
					if ((bc->curr_bytes + add_bytes) >= BUFLEN)
					{
						ret = flush_buffer(FORWARD, bc);
						if (ret)
						{
							logerror("%s:%d Problem writing buffer to file.\n", __func__, __LINE__);
							return 1;
						}
					}
					bc->curr_bytes += add_bytes;
					t = malloc(add_bytes + 1u);
					if (!t)
					{
						logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
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

	/* Deallocate memory */
	free(skip);

	return 0;
}
