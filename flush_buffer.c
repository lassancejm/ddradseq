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
#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "khash.h"
#include "ddradseq.h"

extern int errno;

int flush_buffer(int orient, BARCODE *bc)
{
	char *filename = strdup(bc->outfile);
	char *pch = NULL;
	char *errstr = NULL;
	unsigned char *buffer = (unsigned char*)bc->buffer;
	unsigned char *cbuffer = NULL;
	int ret = 0;
	int fd;
	const int window_bits = 15;
	const int GZIP_ENCODING = 16;
	size_t len = bc->curr_bytes;
	ssize_t w = 0;
	struct flock fl = {F_WRLCK, SEEK_SET, 0, 0, 0};
	struct flock fl2;
	z_stream strm;
	mode_t mode;

	/* Allocate deflate state */
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	ret = deflateInit2(&strm, 4, Z_DEFLATED, window_bits | GZIP_ENCODING,
	                   8, Z_FILTERED);
	if (ret != Z_OK)
		return 1;

	fl.l_pid = getpid();
	memset(&fl2, 0, sizeof(struct flock));

	/* Set permissions if new output file needs to be created */
	mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Convert forward output file name to reverse */
	if (orient == REVERSE)
	{
		pch = strstr(filename, ".R1.fq.gz");
		strncpy(pch, ".R2", 3);
	}

	/* Get output file descriptor */
	fd = open(filename, O_WRONLY | O_CREAT | O_APPEND, mode);
	if (fd < 0)
	{
		errstr = strerror(errno);
		logerror("%s:%d Unable to open output file \'%s\': %s.\n", __func__,
		         __LINE__, filename, errstr);
		return 1;
	}

	/* Test if output file has lock */
	fcntl(fd, F_GETLK, &fl2);
	if (fl2.l_type != F_UNLCK)
	{
		errstr = strerror(errno);
		logerror("%s:%d File \'%s\' is locked for writing by process %zu: %s.\n",
		         __func__, __LINE__, filename, fl.l_pid, errstr);
		return 1;
	}
	else
	{
		/* Set lock on output file */
	    if (fcntl(fd, F_SETLKW, &fl) == -1)
		{
			errstr = strerror(errno);
			logerror("%s:%d Failed to set lock on file \'%s\': %s.\n", __func__,
			         __LINE__, filename, errstr);
			return 1;
		}
	}

	/* Deflate buffer */
	cbuffer = malloc(len);
	if (UNLIKELY(!cbuffer))
	{
		logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 1;
	}

	strm.avail_in = len;
	strm.next_in = buffer;
	strm.avail_out = len;
	strm.next_out = cbuffer;
	ret = deflate(&strm, Z_FINISH);
	if (ret == Z_STREAM_ERROR)
		return 1;

	/* Dump compressed buffer to file */
	w = write(fd, cbuffer, len);
	if (w < 0)
	{
		errstr = strerror(errno);
		logerror("%s:%d Problem writing to output file \'%s\': %s.\n", __func__,
		         __LINE__, filename, errstr);
		return 1;
	}

	/* Deallocate deflate state */
	(void)deflateEnd(&strm);

	/* Reset buffer */
	bc->curr_bytes = 0;
	memset(bc->buffer, 0, BUFLEN);
	bc->buffer[0] = '\0';

	/* Deallocate memory for compressed buffer */
	free(cbuffer);

	/* Unlock output file */
	fl.l_type = F_UNLCK;
	fcntl(fd, F_SETLK, &fl);

	/* Close output file */
	close(fd);

	/* Free allocated memory */
	free(filename);

	return 0;
}
