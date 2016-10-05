/* file: flush_buffer.c
 * description: Dump full buffer to file
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
flush_buffer(int orient, BARCODE *bc)
{
	char *filename = NULL;
	char *buffer = NULL;
	char *pch = NULL;
	int ret = 0;
	size_t len = 0;
	gzFile out;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Get filename */
	filename = strdup(bc->outfile);

	/* Get buffer */
	buffer = bc->buffer;

	/* Get current size of buffer */
	len = bc->curr_bytes;

	/* Convert forward output file name to reverse */
	if (orient == REVERSE)
	{
		pch = strstr(filename, ".R1.fq.gz");
		strncpy(pch, ".R2", 3);
	}

	/* Open output fastQ file stream */
	if ((out = gzopen(filename, "ab")) == NULL)
	{
		fprintf(stderr, "ERROR: Unable to open output file \'%s\'.\n", filename);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Unable to open output file "
		        "\'%s\'.\n", timestr, __func__, __LINE__, filename);
		free(filename);
		return 1;
	}

	/* Dump buffer to file */
	if ((ret = gzwrite(out, buffer, len)) != (int)len)
	{
		fprintf(stderr, "ERROR: Problem writing to output file \'%s\'.\n", filename);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Problem writing to output file "
		        "\'%s\'.\n", timestr, __func__, __LINE__, filename);
		free(filename);
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
