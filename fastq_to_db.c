/* file fastq_to_db.c
 * brief Read in a fastQ database from one or two input files
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

#define MAX_LINE_LENGTH 400
#define BSIZE 4000
#define KEYLEN 31

khash_t(fastq) *
fastq_to_db(char *filename)
{
	char buf[BSIZE][MAX_LINE_LENGTH];
	char *idline = NULL;
	char seps[] = ": ";
	char *tok = NULL;
	char *mkey = NULL;
	char *r = NULL;
	int a = 0;
	int tile = 0;
	int xpos = 0;
	int ypos = 0;
	size_t z = 0;
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

	/* Initialize fastQ hash */
	h = kh_init(fastq);

	/* Open the fastQ input stream */
	if ((in = gzopen(filename, "rb")) == Z_NULL)
	{
		fprintf(stderr, "ERROR: Failed to open input fastQ file %s.\n", filename);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to open input fastQ file %s.\n",
		        timestr, __func__, __LINE__, filename);
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
			if (gzgets(in, buf[lc], MAX_LINE_LENGTH) == Z_NULL) break;
		}

		/* Iterate through lines in the buffer */
		for (l = 0; l < lc; l++)
		{
			/* We are at the end of one fastQ entry */
			if (l % 4 == 3)
			{
				/* Allocate memory for new fastQ entry */
				if ((e = malloc(sizeof(FASTQ))) == NULL)
				{
					fputs("ERROR: Memory allocation failure.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
					        timestr, __func__, __LINE__);
					kh_destroy(fastq, h);
					return NULL;
				}

				/* Parse entry identifier */
				pos = strcspn(buf[l - 3], "\n");
				buf[l - 3][pos] = '\0';
				strl = strlen(&buf[l - 3][1]);
				if ((e->id = malloc(strl + 1u)) == NULL)
				{
					fputs("ERROR: Memory allocation failure.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
					        timestr, __func__, __LINE__);
					free(e);
					kh_destroy(fastq, h);
					return NULL;
				}
				strcpy(e->id, &buf[l - 3][1]);

				/* Construct fastQ hash key */
				if ((idline = malloc(strl + 1u)) == NULL)
				{
					fputs("ERROR: Memory allocation failure.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
					        timestr, __func__, __LINE__);
					free(e->id);
					free(e);
					kh_destroy(fastq, h);
					return NULL;
				}
				strcpy(idline, &buf[l - 3][1]);

				/* Get instrument name */
				if ((tok = strtok_r(idline, seps, &r)) == NULL)
				{
					fputs("ERROR: Parsing ID line failed.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing ID line failed.\n",
					        timestr, __func__, __LINE__);
					free(e->id);
					free(e);
					free(idline);
					kh_destroy(fastq, h);
					return NULL;
				}

				/* Get run number, flow cell ID, and lane number */
				for (z = 0; z < 3; z++)
				{
					if ((tok = strtok_r(NULL, seps, &r)) == NULL)
					{
						fputs("ERROR: Parsing ID line failed.\n", stderr);
						fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing ID line failed.\n",
						        timestr, __func__, __LINE__);
						free(e->id);
						free(e);
						free(idline);
						kh_destroy(fastq, h);
						return NULL;
					}
				}

				/* Get tile number, x, and y coordinate */
				if ((tok = strtok_r(NULL, seps, &r)) == NULL)
				{
					fputs("ERROR: Parsing ID line failed.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing ID line failed.\n",
					        timestr, __func__, __LINE__);
					free(e->id);
					free(e);
					free(idline);
					kh_destroy(fastq, h);
					return NULL;
				}
				tile = atoi(tok);
				if ((tok = strtok_r(NULL, seps, &r)) == NULL)
				{
					fputs("ERROR: Parsing ID line failed.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing ID line failed.\n",
					        timestr, __func__, __LINE__);
					free(e->id);
					free(e);
					free(idline);
					kh_destroy(fastq, h);
					return NULL;
				}
				xpos = atoi(tok);
				if ((tok = strtok_r(NULL, seps, &r)) == NULL)
				{
					fputs("ERROR: Parsing ID line failed.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing ID line failed.\n",
					        timestr, __func__, __LINE__);
					free(e->id);
					free(e);
					free(idline);
					kh_destroy(fastq, h);
					return NULL;
				}
				ypos = atoi(tok);

				/* Construct the hash key */
				if ((mkey = malloc(KEYLEN)) == NULL)
				{
					fputs("ERROR: Memory allocation failure.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
					        timestr, __func__, __LINE__);
					free(e->id);
					free(e);
					free(idline);
					kh_destroy(fastq, h);
					return NULL;
				}
				sprintf(mkey, "%010d%010d%010d", tile, xpos, ypos);
				k = kh_put(fastq, h, mkey, &a);
				if (!a) free(mkey);

				/* Parse DNA sequence */
				pos = strcspn(buf[l - 2], "\n");
				buf[l - 2][pos] = '\0';
				strl = strlen(&buf[l - 2][0]);
				if ((e->seq = malloc(strl + 1u)) == NULL)
				{
					fputs("ERROR: Memory allocation failure.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
					        timestr, __func__, __LINE__);
					free(e->id);
					free(e);
					free(idline);
					kh_destroy(fastq, h);
					return NULL;
				}
				strcpy(e->seq, &buf[l - 2][0]);

				/* Parse quality sequence */
				pos = strcspn(buf[l], "\n");
				buf[l][pos] = '\0';
				strl = strlen(&buf[l][0]);
				if ((e->qual = malloc(strl + 1u)) == NULL)
				{
					fputs("ERROR: Memory allocation failure.\n", stderr);
					fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
					        timestr, __func__, __LINE__);
					free(e->id);
					free(e);
					free(idline);
					free(e->seq);
					kh_destroy(fastq, h);
					return NULL;
				}
				strcpy(e->qual, &buf[l][0]);

				/* Add to database */
				kh_value(h, k) = e;

				/* Free allocated memory */
				free(idline);
			}
		}

		/* If we are at the end of the file */
		if (lc < BSIZE) break;
	}

	/* Close input stream */
	gzclose(in);

	return h;
}
