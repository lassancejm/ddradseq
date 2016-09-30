/* file: ddradseq.c
 * description: Entry point for the ddradseq program
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "ddradseq.h"

/* Function prototypes */
extern int parse_main(CMD*);
extern int trimend_main(CMD*);
extern int pair_main(CMD*);
void usage(void);

int
main(int argc, char *argv[])
{
	char logfname[] = "/ddradseq.log\0";
	int ret = 0;
	CMD *cp = NULL;

	/* Get current working directory */
	if (getcwd(logfile, sizeof(logfile)) == NULL)
	{
		fputs("ERROR: Failed to get current working directory.\n", stderr);
		return 1;
	}
	strcat(logfile, logfname);

	/* Open log file for writing */
	lf = fopen(logfile, "a");

	/* Parse the command line options */
	if ((cp = parse_cmdline(argc, argv)) == NULL)
	{
		usage();
		return 1;
	}

	/* Initialize the log files */
	log_init(ret);

	/* Run the pipeline stages */
	if (strcmp(cp->mode, "parse") == 0 || strcmp(cp->mode, "all") == 0)
		if ((ret = parse_main(cp)) != 0)
			return 1;
	if (strcmp(cp->mode, "pair") == 0 || strcmp(cp->mode, "all") == 0)
		if ((ret = pair_main(cp)) != 0)
			return 1;
	if (strcmp(cp->mode, "trimend") == 0 || strcmp(cp->mode, "all") == 0)
		if ((ret = pair_main(cp)) != 0)
			return 1;

	/* Close logfile output stream */
	fclose(lf);

	free_cmdline(cp);

	return 0;
}

void
usage(void)
{
	fputs("Usage: ddradseq [OPTIONS] [INPUT DIRECTORY]\n\n", stderr);
	fputs("Parse fastQ file into separate files by flowcell, barcode and/or index\n\n", stderr);
	fputs("Mandatory arguments to long options are mandatory for short options too.\n", stderr);
	fputs(" -m  --mode=STR       Run mode of ddradseq program\n", stderr);
	fputs("                      Valid modes are:\n", stderr);
	fputs("	                     parse     Parses input fastQ by standard Illumina and custom adapter\n", stderr);
	fputs("	                     pair      Aligns mated pairs in two fastQ input files\n", stderr);
	fputs("	                     trimend   Trims the 3\' end of fastQ reverse sequences\n", stderr);
	fputs(" -c, --csv=FILE       CSV file with index and barcode labels\n", stderr);
	fputs(" -o, --out=DIR        Parent directory to write output files                 [default: same as input fastQ]\n", stderr);
	fputs(" -d, --dist           Edit distance for barcode matching                     [default: 1]\n", stderr);
	fputs(" -s, --score          Alignment score to consider sequences properly aligned [default: 100]\n", stderr);
	fputs(" -h, --help           Display this help message\n", stderr);
	fputs(" -v, --version        Print software version number and exit\n\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
