/* file: parse_fastq.c
 * description: Main code to parse a fastQ file according to index sequence
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
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
parse_fastq (char *filename, khash_t(pool_hash) *h)
{
	char *r = NULL;
	char *p = NULL;
	char buffer[BUFLEN];
	size_t numlines = 0;
	size_t bytes_read = 0;
	size_t buff_rem = 0;
	gzFile fin;

	/* Open input file */
	fin = gzopen(filename, "rb");
	if (fin == NULL)
	{
		fprintf(stderr, "Error: unable to open file \'%s\'.\n", filename);
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
		p = &buffer[0];
		numlines = count_lines(p);
		r = clean_buffer(p, &numlines);
		parse_buffer(p, numlines, h);
		buff_rem = reset_buffer(p, r);

		/* Check if we are at the end of file */
		if (gzeof(fin)) break;
	}

	/* Close input file */
	gzclose(fin);

	return 0;
}
