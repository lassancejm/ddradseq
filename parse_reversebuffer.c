/* file: parse_reversebuffer.c
 * description: Parses reverse fastQ entries in the buffer
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

int parse_reversebuffer(char *buff, const size_t nl, const khash_t(pool_hash) *h,
                        const khash_t(mates) *m)
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
	unsigned char *skip = NULL;
	int ret = 0;
	int x = 0;
	int tile = 0;
	int xpos = 0;
	int ypos = 0;
	size_t add_bytes = 0;
	size_t l = 0;
	size_t ll = 0;
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
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
						logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					idline = strdup(q);
					if (!idline)
					{
						logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}

					/* Parse Illumina identifier line */
					for (tok = strtok_r(copy, seps, &r), x = 0; tok != NULL; tok = strtok_r(NULL, seps, &r), x++)
					{
						switch (x)
						{
							case 2:
								flowcell_ID = strdup(tok);
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

					/* Check parsing results */
					if (!flowcell_ID)
					{
						logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					if (!index_sequence)
					{
						logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}

					/* Lookup flow cell identifier */
					i = kh_get(pool_hash, h, flowcell_ID);

					/* Flow cell identifier is not present in database */
					if (i == kh_end(h))
					{
						fprintf(lf, "[ddradseq: %s] WARNING -- Hash lookup failure using key %s.\n", timestr, flowcell_ID);
						fprintf(lf, "[ddradseq: %s] WARNING -- Skipping sequence: %s\n", timestr, idline);
						skip[l+1] = 1;
						skip[l+2] = 1;
						skip[l+3] = 1;
						free(idline);
						free(copy);
						free(flowcell_ID);
						break;
					}
					else
						p = kh_value(h, i);

					/* Lookup pool identifier */
					j = kh_get(pool, p, index_sequence);
					pl = kh_value(p, j);
					b = pl->b;

					/* Retrieve barcode sequence of mate */
					mkey = malloc(KEYLEN);
					if (!mkey)
					{
						logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					sprintf(mkey, "%010d%010d%010d", tile, xpos, ypos);
					mk = kh_get(mates, m, mkey);
					if (mk == kh_end(m))
					{
						fprintf(lf, "[ddradseq: %s] WARNING -- Hash lookup failure using key %s.\n", timestr, mkey);
						fprintf(lf, "[ddradseq: %s] WARNING -- Skipping sequence: %s\n", timestr, idline);
						skip[l+1] = 1;
						skip[l+2] = 1;
						skip[l+3] = 1;
						free(mkey);
						free(idline);
						free(copy);
						free(flowcell_ID);
						break;
					}
					barcode_sequence = kh_value(m, mk);
					free(mkey);

					/* Get the barcode entry of read's mate */
					k = kh_get(barcode, b, barcode_sequence);
					if (k != kh_end(b))
						bc = kh_value(b, k);
					else
					{
						skip[l+1] = 1;
						skip[l+2] = 1;
						skip[l+3] = 1;
						free(idline);
					}

					/* Free memory */
					free(flowcell_ID);
					free(index_sequence);
					free(copy);
					break;
				case 1:
					/* Sequence line */
					dna_sequence = strdup(q);
					if (!dna_sequence)
					{
						logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					break;
				case 2:
					/* Quality identifier line */
					break;
				case 3:
					/* Quality sequence line */
					qual_sequence = strdup(q);
					if (!qual_sequence)
					{
						logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					add_bytes = strlen(idline) + strlen(dna_sequence) +
								strlen(qual_sequence) + 5u;
					char *t = NULL;
					if ((bc->curr_bytes + add_bytes) >= BUFLEN)
					{
						ret = flush_buffer(REVERSE, bc);
						if (ret)
						{
							logerror("%s:%d Problem writing to file.\n", __func__, __LINE__);
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
					sprintf(t, "%s\n%s\n+\n%s\n", idline, dna_sequence, qual_sequence);
					strcat(bc->buffer, t);

					/* Free alloc'd memory for fastQ entry */
					free(t);
					free(idline);
					free(dna_sequence);
					free(qual_sequence);
					break;
			}
		}
		q += ll + 1u;
	}

	/* Deallocate memory */
	free(skip);

	return 0;
}
