/* file pair_mates.c
 * description: Pair mates in two fastQ files
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

int pair_mates(const char *filename, const khash_t(fastq) *h, const char *ffor,
               const char *frev)
{
	char **buf = NULL;
	char *idline = NULL;
	char seps[] = ": ";
	char *tok = NULL;
	char *mkey = NULL;
	char *r = NULL;
	int tile = 0;
	int i = 0;
	int xpos = 0;
	int ypos = 0;
	size_t x = 0;
	size_t l = 0;
	size_t lc = 0;
	size_t pos = 0;
	size_t strl = 0;
	khint_t k = 0;
	gzFile in;
	gzFile fout;
	gzFile rout;
	FASTQ *e = NULL;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Allocate memory for buffer from heap */
	buf = malloc(BSIZE * sizeof(char*));
	if (UNLIKELY(!buf))
	{
		logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 1;
	}
	for (i = 0; i < BSIZE; i++)
	{
		buf[i] = malloc(MAX_LINE_LENGTH);
		if (UNLIKELY(!buf[i]))
		{
			logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
	}

	/* Open the fastQ input stream */
	in = gzopen(filename, "rb");
	if (!in)
	{
		logerror("%s:%d Failed to open input fastQ file \'%s\'.\n", __func__,
		         __LINE__, filename);
		return 1;
	}

	/* Open the output fastQ file streams */
	fout = gzopen(ffor, "wb");
	if (!fout)
	{
		logerror("%s:%d Failed to open forward output fastQ file \'%s\'.\n",
		         __func__, __LINE__, ffor);
		gzclose(in);
		return 1;
	}

	rout = gzopen(frev, "wb");
	if (!rout)
	{
		logerror("%s:%d Failed to open reverse output fastQ file \'%s\'.\n",
		         __func__, __LINE__, frev);
		gzclose(in);
		gzclose(fout);
		return 1;
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
				/* Parse entry identifier */
				pos = strcspn(buf[l-3], "\n");
				buf[l-3][pos] = '\0';
				strl = strlen(&buf[l-3][1]);

				/* Construct fastQ hash key */
				idline = malloc(strl + 1u);
				if (UNLIKELY(!idline))
				{
					logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return 1;
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

				/* Construct hash key */
				mkey = malloc(KEYLEN);
				if (UNLIKELY(!mkey))
				{
					logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return 1;
				}
				sprintf(mkey, "%010d%010d%010d", tile, xpos, ypos);
				k = kh_get(fastq, h, mkey);
				if (k != kh_end(h))
					e = kh_value(h, k);
				free(mkey);
				free(idline);

				if (e != NULL)
				{
					/* Parse DNA sequence */
					pos = strcspn(buf[l-2], "\n");
					buf[l-2][pos] = '\0';

					/* Parse quality sequence */
					pos = strcspn(buf[l], "\n");
					buf[l][pos] = '\0';

					/* Need to construct output file streams */
					gzprintf(fout, "@%s\n%s\n+\n%s\n", e->id, e->seq, e->qual);
					gzprintf(rout, "@%s\n%s\n+\n%s\n", &buf[l-3][1], &buf[l-2][0], &buf[l][0]);
				}
			}
		}

		/* If we are at the end of the file */
		if (lc < BSIZE)
			break;
	}

	/* Free memory for buffer to heap */
	for (i = 0; i < BSIZE; i++)
		free(buf[i]);
	free(buf);

	/* Close input stream */
	gzclose(in);
	gzclose(fout);
	gzclose(rout);

	return 0;
}
