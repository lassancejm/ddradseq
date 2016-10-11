/* file fastq_to_db.c
 * brief Populates a fastQ database from fastQ input file
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "ddradseq.h"
#include "khash.h"

khash_t(fastq) *fastq_to_db(const char *filename)
{
	char **buf = NULL;
	char *idline = NULL;
	char seps[] = ": ";
	char *tok = NULL;
	char *mkey = NULL;
	char *r = NULL;
	int a = 0;
	int i = 0;
	int tile = 0;
	int xpos = 0;
	int ypos = 0;
	size_t x = 0;
	size_t l = 0;
	size_t lc = 0;
	size_t pos = 0;
	size_t strl = 0;
	khint_t k = 0;
	gzFile in;
	FASTQ *e = NULL;
	khash_t(fastq) *h = NULL;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Allocate memory for buffer from heap */
	buf = malloc(BSIZE * sizeof(char*));
	if (UNLIKELY(!buf))
	{
		logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return NULL;
	}
	for (i = 0; i < BSIZE; i++)
	{
		buf[i] = malloc(MAX_LINE_LENGTH);
		if (UNLIKELY(!buf[i]))
		{
			logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return NULL;
		}
	}

	/* Initialize fastQ hash */
	h = kh_init(fastq);

	/* Open the fastQ input stream */
	in = gzopen(filename, "rb");
	if (!in)
	{
		logerror("%s:%d Failed to open input fastQ file %s.\n", __func__,
			     __LINE__, filename);
		return NULL;
	}

	/* Enter data from the fastQ input file into the database */
	while (1)
	{
		/* Fill up the buffer */
		for (lc = 0; lc < BSIZE; lc++)
		{
			/* Clear the buffer */
			memset(buf[lc], 0, MAX_LINE_LENGTH);

			/* Get line from the fastQ input stream */
			if (gzgets(in, buf[lc], MAX_LINE_LENGTH) == Z_NULL)
				break;
		}

		/* Iterate through lines in the buffer */
		for (l = 0; l < lc; l++)
		{
			/* We are at the end of one fastQ entry */
			if (l % 4 == 3)
			{
				/* Allocate memory for new fastQ entry */
				e = malloc(sizeof(FASTQ));
				if (UNLIKELY(!e))
				{
					logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return NULL;
				}

				/* Parse entry identifier */
				pos = strcspn(buf[l-3], "\n");
				buf[l-3][pos] = '\0';
				strl = strlen(&buf[l-3][1]);
				e->id = malloc(strl + 1u);
				if (UNLIKELY(!e->id))
				{
					logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return NULL;
				}
				strcpy(e->id, &buf[l-3][1]);

				/* Construct fastQ hash key */
				idline = malloc(strl + 1u);
				if (UNLIKELY(!idline))
				{
					logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return NULL;
				}
				strcpy(idline, &buf[l-3][1]);

				/* Parse Illumina identifier line */
				for (tok = strtok_r(idline, seps, &r), x = 0; tok != NULL; tok = strtok_r(NULL, seps, &r), x++)
				{
					switch (x)
					{
						case 4:
							tile = atoi(tok);
							break;
						case 5:
							xpos = atoi(tok);
							break;
						case 6:
							ypos = atoi(tok);
							break;
					}
				}

				/* Construct the hash key */
				mkey = malloc(KEYLEN);
				if (UNLIKELY(!mkey))
				{
					logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return NULL;
				}
				sprintf(mkey, "%010d%010d%010d", tile, xpos, ypos);
				k = kh_put(fastq, h, mkey, &a);
				if (!a)
					free(mkey);

				/* Parse DNA sequence */
				pos = strcspn(buf[l-2], "\n");
				buf[l-2][pos] = '\0';
				strl = strlen(&buf[l-2][0]);
				e->seq = malloc(strl + 1u);
				if (UNLIKELY(!e->seq))
				{
					logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return NULL;
				}
				strcpy(e->seq, &buf[l-2][0]);

				/* Parse quality sequence */
				pos = strcspn(buf[l], "\n");
				buf[l][pos] = '\0';
				strl = strlen(&buf[l][0]);
				e->qual = malloc(strl + 1u);
				if (UNLIKELY(!e->qual))
				{
					logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return NULL;
				}
				strcpy(e->qual, &buf[l][0]);

				/* Add to database */
				kh_value(h, k) = e;

				/* Free allocated memory from heap */
				free(idline);
			}
		}

		/* If we are at the end of the file */
		if (lc < BSIZE) break;
	}

	/* Free memory for buffer to heap */
	for (i = 0; i < BSIZE; i++)
		free(buf[i]);
	free(buf);

	/* Close input stream */
	gzclose(in);

	return h;
}
