/* file: flush_buffer.c
 * description: Dump full buffer to file
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "khash.h"
#include "ddradseq.h"

int flush_buffer(BARCODE *bc)
{
	char *filename = bc->outfile;
	char *buffer = bc->buffer;
	int ret = 0;
	size_t len = bc->curr_bytes;
	gzFile out;

	/* Open input database text file stream */
	if ((out = gzopen(filename, "ab")) == NULL)
	{
		fprintf(stderr, "Error opening output file \'%s\'.\n", filename);
		return 1;
	}

	/* Dump buffer to file */
	ret = gzwrite(out, buffer, len);
	if (ret != (int)len)
	{
		fprintf(stderr, "Error writing to output file \'%s\'.\n", filename);
		return 1;
	}

	/* Reset buffer */
	bc->curr_bytes = 0;
	memset(bc->buffer, 0, BUFLEN);
	bc->buffer[0] = '\0';

	/* Close output file */
	gzclose(out);

	return 0;
}

int flush_pbuffer(POOL *pl)
{
	char *filename = pl->poutfile;
	char *buffer = pl->pbuffer;
	int ret = 0;
	size_t len = pl->pcurr_bytes;
	gzFile out;

	/* Open input database text file stream */
	if ((out = gzopen(filename, "ab")) == NULL)
	{
		fprintf(stderr, "Error opening output file \'%s\'.\n", filename);
		return 1;
	}

	/* Dump buffer to file */
	ret = gzwrite(out, buffer, len);
	if (ret != (int)len)
	{
		fprintf(stderr, "Error writing to output file \'%s\'.\n", filename);
		return 1;
	}

	/* Reset buffer */
	pl->pcurr_bytes = 0;
	memset(pl->pbuffer, 0, BUFLEN);
	pl->pbuffer[0] = '\0';

	/* Close output file */
	gzclose(out);

	return 0;
}

