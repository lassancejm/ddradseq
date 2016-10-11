/* file: flush_buffer.c
 * description: Dumps a full buffer to file
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#include "khash.h"
#include "ddradseq.h"

extern int errno;

int flush_buffer(int orient, BARCODE *bc)
{
	char *filename = NULL;
	char *buffer = NULL;
	char *pch = NULL;
	char *errstr = NULL;
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
	out = gzopen(filename, "ab");
	if (!out)
	{
		errstr = strerror(errno);
		logerror("%s:%d Unable to open output file \'%s\': %s.\n", __func__,
		         __LINE__, filename, errstr);
		return 1;
	}

	/* Dump buffer to file */
	ret = gzwrite(out, buffer, len);
	if (ret == 0)
	{
		logerror("%s:%d Problem writing to output file \'%s\'.\n", __func__,
		         __LINE__, filename);
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
