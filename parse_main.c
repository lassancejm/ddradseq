/* file: parse_main.c
 * description: Entry point for the parse modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include "ddradseq.h"

void parse_usage(void);

int
parse_main(int argc, char *argv[])
{
	int ret = 0;
	CMD *cp = NULL;
	khash_t(pool_hash) *h = NULL;
	khash_t(mates) *m = NULL;

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "parse")) == NULL)
	{
		parse_usage();
		ret = 1;
	}

	/* Read CSV database into memory */
    if (ret == 0)
        h = read_csv(cp);

	/* Check for write permissions on parent of output directory */
	if ((ret == 0) && (h != NULL))
		ret = check_directories(cp, h);

	/* Initialize hash for mate pair information */
	m = kh_init(mates);

	/* Read the fastQ input files */
	if ((ret == 0) && (h != NULL))
	{
		parse_fastq(FORWARD, cp->forfastq, h, m, cp->trim_barcode);
		parse_fastq(REVERSE, cp->revfastq, h, m, cp->trim_barcode);
	}

	/* Deallocate memory */
    free_db(h);
    kh_destroy(mates, m);
	if (cp)
		free_cmdline(cp);

	return ret;
}

void
parse_usage(void)
{
	fputs("Usage : ddradseq parse [OPTIONS] [FASTQ.R1] [FASTQ.R2]\n\n", stderr);
	fputs("Parse fastQ file into separate files by flowcell, barcode and/or index\n", stderr);
	fputs("If a custom barcode is used it is automatically trimmed from the 5\' end of the forward sequence.\n\n", stderr);
	fputs("Available options\n", stderr);
	fputs("  -i  FILE   CSV file with index and barcode labels\n", stderr);
	fputs("  -o  DIR    Parent directory to write output files [default: same as input fastQ]\n", stderr);
	fputs("  -n  INT    Number of threads for concurrency [default: 1]\n", stderr);
	fputs("  -t         Trim barcodes from sequences\n", stderr);
	fputs("  -h         Display this help message\n\n", stderr);
}
