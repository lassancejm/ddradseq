/* file: pair_main.c
 * description: Entry point for the pair modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include "ddradseq.h"

void pair_usage(void);

int
pair_main(int argc, char *argv[])
{
	int ret = 0;
	CMD *cp = NULL;

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "pair")) == NULL)
	{
		pair_usage();
		ret = 1;
	}

	if (ret == 0)
		pair_usage();

	/* Deallocate memory */
	if (cp)
		free_cmdline(cp);

	return ret;
}

void
pair_usage(void)
{
	fputs("Usage : ddradseq pair [OPTIONS] [FASTQ.R1] [FASTQ.R2]\n\n", stderr);
	fputs("Aligns mated pairs in two fastQ files.\n\n", stderr);
	fputs("Available options\n", stderr);
	fputs("  -o  DIR    Parent directory to write output files [default: same as input fastQ]\n", stderr);
	fputs("  -t  INT    Number of threads for concurrency [default: 1]\n", stderr);
	fputs("  -h         Display this help message\n\n", stderr);
}

