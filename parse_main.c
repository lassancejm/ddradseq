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
void deallocate_mem(CMD*, khash_t(pool_hash)*, khash_t(mates)*);

int
parse_main(int argc, char *argv[])
{
	int ret = EXIT_SUCCESS;
	CMD *cp = NULL;
	khash_t(pool_hash) *h = NULL;
	khash_t(mates) *m = NULL;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv, "parse")) == NULL)
	{
		parse_usage();
		return EXIT_FAILURE;
	}

	/* Read CSV database into memory */
	h = read_csv(cp);
	if (h == NULL)
	{
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to read CSV database "
		        "into memory.\n", timestr, __func__, __LINE__);
		deallocate_mem(cp, NULL, NULL);
		return EXIT_FAILURE;
	}

	/* Check for write permissions on parent of output directory */
	if ((ret = check_directories(cp, h)) == EXIT_FAILURE)
	{
		deallocate_mem(cp, h, NULL);
		return EXIT_FAILURE;
	}

	/* Initialize hash for mate pair information */
	if ((m = kh_init(mates)) == NULL)
	{
		deallocate_mem(cp, h, NULL);
		return EXIT_FAILURE;
	}

	/* Read the forward fastQ input file */
	if ((ret = parse_fastq(FORWARD, cp->forfastq, h, m, cp->dist)) == EXIT_FAILURE)
	{
		deallocate_mem(cp, h, m);
		return EXIT_FAILURE;
	}

	/* Read the reverse fastQ input file */
	if ((ret = parse_fastq(REVERSE, cp->revfastq, h, m, cp->dist)) == EXIT_FAILURE)
	{
		deallocate_mem(cp, h, m);
		return EXIT_FAILURE;
	}

	/* Deallocate memory */
	deallocate_mem(cp, h, m);

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log */
	fprintf(lf, "[ddradseq: %s] INFO -- Parse step of pipeline is complete.\n", timestr);

	return EXIT_SUCCESS;
}

void
deallocate_mem(CMD *cp, khash_t(pool_hash) *h, khash_t(mates) *m)
{
	if (h) free_db(h);
	if (m) kh_destroy(mates, m);
	if (cp) free_cmdline(cp);	
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
