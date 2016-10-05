/* file: parse_fastq.c
 * description: Main code to parse a fastQ file by index sequence
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "khash.h"
#include "ddradseq.h"

int
parse_fastq(int orient, char *filename, khash_t(pool_hash) *h, khash_t(mates) *m, int dist)
{
	char *r = NULL;
	char *q = NULL;
	char buffer[BUFLEN];
	int ret = 0;
	size_t numlines = 0;
	size_t bytes_read = 0;
	size_t buff_rem = 0;
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khash_t(barcode) *b = NULL;
	khash_t(pool) *p = NULL;
	BARCODE *bc = NULL;
	POOL *pl = NULL;
	gzFile fin;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log */
	fprintf(lf, "[ddradseq: %s] INFO -- Parsing fastQ file \'%s\'.\n", timestr,
	        filename);

	/* Open input file */
	fin = gzopen(filename, "rb");
	if (fin == NULL)
	{
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Unable to open file \'%s\'.\n",
		        timestr, __func__, __LINE__, filename);
		return 1;
	}

	/* Initialize buffer */
	memset(buffer, 0, sizeof(buffer));
	r = &buffer[0];

	/* Iterate through blocks from input fastQ file */
	while (1)
	{
		/* Read block from file into input buffer */
		bytes_read = gzread(fin, &buffer[buff_rem], BUFLEN - buff_rem - 1);

		/* Set null terminating character on input buffer */
		buffer[bytes_read + buff_rem] = '\0';

		/* Count lines in input buffer */
		q = &buffer[0];
		numlines = count_lines(q);
		r = clean_buffer(q, &numlines);
		if (orient == FORWARD)
		{
			if ((ret = parse_forwardbuffer(q, numlines, h, m, dist)) != 0)
				return 1;
		}
		else
		{
			if ((ret = parse_reversebuffer(q, numlines, h, m)) != 0)
				return 1;
		}
		buff_rem = reset_buffer(q, r);

		/* Check if we are at the end of file */
		if (gzeof(fin)) break;
	}

	/* Flush remaining data in buffers */
	for (i = kh_begin(h); i != kh_end(h); i++)
	{
		if (kh_exist(h, i))
		{
			p = kh_value(h, i);
			for (j = kh_begin(p); j != kh_end(p); j++)
			{
				if (kh_exist(p, j))
				{
					pl = kh_value(p, j);
					b = pl->b;
					for (k = kh_begin(b); k != kh_end(b); k++)
					{
						if (kh_exist(b, k))
						{
							bc = kh_value(b, k);
							if (bc->curr_bytes > 0)
							{
								if ((ret = flush_buffer(orient, bc)) != 0)
								{
									fputs("ERROR: Problem writing buffer to file.\n", stderr);
									fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Problem writing buffer to file.\n",
									        timestr, __func__, __LINE__);
									return 1;
								}
							}
						}
					}
				}
			}
		}
	}

	/* Close input file */
	gzclose(fin);

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log */
	fprintf(lf, "[ddradseq: %s] INFO -- Successfully parsed fastQ file \'%s\'.\n",
	        timestr, filename);

	return 0;
}
