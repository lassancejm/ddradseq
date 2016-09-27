/* file: trimend_main.c
 * description: Entry point for the trimend function
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include "ddradseq.h"

void trimend_usage (void);

int
trimend_main(int argc, char *argv[])
{
	int ret = 0;
	CMD *cp = NULL;

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "trimend")) == NULL)
	{
		trimend_usage();
		ret = 1;
	}

	if (ret == 0)
		trimend_usage();

	/* Deallocate memory */
	if (cp)
		free_cmdline(cp);

	return ret;
}

void
trimend_usage(void)
{
	fputs("Usage : ddradseq trimend [OPTIONS] [FASTQ.R1] [FASTQ.R2]\n\n", stderr);
	fputs("Trims the 5\' end of reverse reads.\n", stderr);
	fputs("Mated pairs are aligned and any overhand is trimmed.\n\n", stderr);
	fputs("Mandatory arguments to long options are mandatory for short options too.\n", stderr);
	fputs(" -o, --out=DIR        Parent directory to write output files    [default: same as input fastQ]\n", stderr);
	fputs(" -n, --threads=INT    Number of threads for concurrency         [default: 1]\n", stderr);
	fputs(" -h, --help           Display this help message\n\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
