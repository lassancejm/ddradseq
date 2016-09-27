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
#include "khash.h"

void pair_usage(void);

int
pair_main(int argc, char *argv[])
{
	int ret = 0;
	CMD *cp = NULL;
	khash_t(fastq) *h = NULL;

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "pair")) == NULL)
	{
		pair_usage();
		ret = 1;
	}

	if (ret == 0)
		pair_usage();

	/* Read forward fastQ file into hash table */
	h = fastq_to_db(cp->forfastq);

	/* Align mated pairs and write to output file*/
	pair_mates(cp->revfastq, h);

	/* Deallocate memory */
	free_pairdb(h);
	if (cp)
		free_cmdline(cp);

	return ret;
}

void
pair_usage(void)
{
	fputs("Usage : ddradseq pair [OPTIONS] [FASTQ.R1] [FASTQ.R2]\n\n", stderr);
	fputs("Aligns mated pairs in two fastQ files.\n\n", stderr);
	fputs("Mandatory arguments to long options are mandatory for short options too.\n", stderr);
	fputs(" -o, --out=DIR        Parent directory to write output files    [default: same as input fastQ]\n", stderr);
	fputs(" -n, --threads=INT    Number of threads for concurrency         [default: 1]\n", stderr);
	fputs(" -h, --help           Display this help message\n\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
