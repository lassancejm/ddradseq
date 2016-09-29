/* file: parse_main.c
 * description: Entry point for the parse modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include "khash.h"
#include "ddradseq.h"

/* Function prototypes */
void parse_usage(void);

int
parse_main(int argc, char *argv[])
{
	CMD *cp = NULL;
	khash_t(pool_hash) *h = NULL;
	khash_t(mates) *m = NULL;

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "parse")) == NULL)
	{
		parse_usage();
		exit(EXIT_FAILURE);
	}

	/* Read CSV database into memory */
	h = read_csv(cp);

	/* Check for write permissions on parent of output directory */
	check_directories(cp, h);

	/* Initialize hash for mate pair information */
	m = kh_init(mates);

	/* Read the fastQ input files */
	parse_fastq(FORWARD, cp->forfastq, h, m, cp->dist);
	parse_fastq(REVERSE, cp->revfastq, h, m, cp->dist);

	/* Deallocate memory */
	free_db(h);
	kh_destroy(mates, m);
	if (cp)
		free_cmdline(cp);

	return 0;
}

void
parse_usage(void)
{
	fputs("Usage : ddradseq parse [OPTIONS] [FASTQ.R1] [FASTQ.R2]\n\n", stderr);
	fputs("Parse fastQ file into separate files by flowcell, barcode and/or index\n\n", stderr);
	fputs("Mandatory arguments to long options are mandatory for short options too.\n", stderr);
	fputs(" -c, --csv=FILE       CSV file with index and barcode labels\n", stderr);
	fputs(" -o, --out=DIR        Parent directory to write output files    [default: same as input fastQ]\n", stderr);
	fputs(" -d, --dist           Edit distance for barcode matching        [default: 1]\n", stderr);
	fputs(" -h, --help           Display this help message\n\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
