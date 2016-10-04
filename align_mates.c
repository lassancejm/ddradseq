/* file align_mates.c
 * description: Align mates in two fastQ files and trim 3' end of reverse sequences
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "ddradseq.h"

/* Define some lengths*/
#define MAX_LINE_LENGTH 400
#define BSIZE 4000

/*Globally scoped variables */
unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,	 4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,	 4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,	 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
char alpha[5] = "ACGTN";

int
align_mates(CMD *cp, char *forin, char *revin, char *forout, char *revout)
{
	char fbuf[BSIZE][MAX_LINE_LENGTH];
	char rbuf[BSIZE][MAX_LINE_LENGTH];
	char mat[25];
	int i = 0;
	int j = 0;
	int k = 0;
	int xtra = KSW_XSTART;
	int sa = 1;
	int sb = 3;
	int gap_open = cp->gapo;
	int gap_extend = cp->gape;
	int min_score = cp->score;
	int count = 0;
	size_t l = 0;
	size_t lc = 0;
	gzFile fin;
	gzFile rin;
	gzFile fout;
	gzFile rout;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Open input forward fastQ file stream */
	if ((fin = gzopen(forin, "rb")) == NULL)
	{
		fprintf(stderr, "ERROR: Failed to open input forward fastQ file \'%s\'.\n", forin);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to open input forward fastQ file \'%s\'.\n",
		        timestr, __func__, __LINE__, forin);
		return 1;
	}

	/* Open input reverse fastQ file stream */
	if ((rin = gzopen(revin, "rb")) == NULL)
	{
		fprintf(stderr, "ERROR: Failed to open input reverse fastQ file \'%s\'.\n", revin);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to open input reverse fastQ file \'%s\'.\n",
		        timestr, __func__, __LINE__, revin);
		gzclose(fin);
		return 1;
	}

	/* Open output forward fastQ file stream */
	if ((fout = gzopen(forout, "wb")) == NULL)
	{
		fprintf(stderr, "ERROR: Failed to open forward output fastQ file \'%s\'.\n", forout);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to open forward output fastQ file \'%s\'.\n",
		        timestr, __func__, __LINE__, forout);
		gzclose(fin);
		gzclose(rin);
		return 1;
	}

	/* Open output reverse fastQ file stream */
	if ((rout = gzopen(revout, "wb")) == NULL)
	{
		fprintf(stderr, "ERROR: Failed to open reverse output fastQ file \'%s\'.\n", revout);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to open reverse output fastQ file \'%s\'.\n",
		        timestr, __func__, __LINE__, revout);
		gzclose(fin);
		gzclose(rin);
		gzclose(fout);
		return 1;
	}

	/* Initialize the scoring matrix */
	for (i = k = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
			mat[k++] = (i == j) ? sa : -sb;

		/* Ambiguous base */
		mat[k++] = 0;
	}

	for (j = 0; j < 5; j++)
		mat[k++] = 0;

	while (1)
	{
		/* Fill up the forward buffer */
		for (lc = 0; lc < BSIZE; lc++)
		{
			/* Clear the buffer */
			memset(fbuf[lc], 0, MAX_LINE_LENGTH);

			/* Get line from the fastQ input stream */
			if (gzgets(fin, fbuf[lc], MAX_LINE_LENGTH) == Z_NULL) break;
		}

		/* Fill up the reverse buffer */
		for (lc = 0; lc < BSIZE; lc++)
		{
			/* Clear the buffer */
			memset(rbuf[lc], 0, MAX_LINE_LENGTH);

			/* Get line from the fastQ input stream */
			if (gzgets(rin, rbuf[lc], MAX_LINE_LENGTH) == Z_NULL) break;
		}

		/* Iterate through lines in the buffers */
		for (l = 0; l < lc; l++)
		{
			/* We are at the end of one fastQ entry */
			if (l % 4 == 3)
			{
				ALIGN_RESULT r;
				char *target;
				char *query;
				target = strdup(&fbuf[l - 2][0]);
				query = revcom(&rbuf[l - 2][0]);
				int tlen = (int)strlen(target);
				target[tlen] = '\0';
				int qlen = (int)strlen(query);

				/* Transform sequences */
				for (i = 0; i < qlen; i++)
					query[i] = seq_nt4_table[(int)query[i]];

				for (i = 0; i < tlen; i++)
					target[i] = seq_nt4_table[(int)target[i]];

				/* Do the alignment */
				r = local_align ((int)qlen,
								 (unsigned char *)query,
								 (int)tlen,
								 (unsigned char *)target,
								 5, mat, gap_open,
								 gap_extend, xtra);
				free(target);
				free(query);

				/* Actually trim the sequence */
				if (r.score >= min_score)
				{
					/* Test trimming criterion */
					if ((r.target_begin == 0) &&
						(r.query_begin > 0))
					{
						int new_end_pos = qlen - r.query_begin;
						char *seq = &rbuf[l - 2][0];
						char *qual = &rbuf[l][0];
						seq[new_end_pos] = '\0';
						qual[new_end_pos] = '\0';
						count++;
					}
				}

				/* Write sequences to file */
				gzprintf(fout, "%s\n%s\n+\n%s\n", &fbuf[l - 3][0], &fbuf[l - 2][0], &fbuf[l][0]);
				gzprintf(rout, "%s\n%s\n+\n%s\n", &rbuf[l - 3][0], &rbuf[l - 2][0], &rbuf[l][0]);
			}
		}

		/* If we are at the end of the file */
		if (lc < BSIZE) break;
	}

	/* Print informational message to logfile */
	get_timestr(&timestr[0]);
	fprintf(lf, "[ddradseq: %s] INFO -- %d sequences trimmed.\n", timestr, count);

	/* Close all file streams */
	gzclose(fin);
	gzclose(rin);
	gzclose(fout);
	gzclose(rout);

	return 0;
}
