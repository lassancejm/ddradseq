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

int flush_buffer(int orient, BARCODE *bc)
{
	char *filename = strdup(bc->outfile);
	char *buffer = bc->buffer;
	char *pch = NULL;
	int ret = 0;
	size_t len = bc->curr_bytes;
	gzFile out;

	/* Convert forward output file name to reverse */
	if (orient == REVERSE)
	{
		pch = strstr(filename, ".R1.fq.gz");
		strncpy(pch, ".R2", 3);
	}

	/* Open output fastQ file stream */
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

	/* Free allocated memory */
	free(filename);

	return 0;
}
