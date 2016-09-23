/* file: ddradseq.c
 * description: Entry point for the ddradseq pipeline worker program
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

extern int parse_main(int argc, char *argv[]);
extern int trimend_main(int argc, char *argv[]);
extern int pair_main(int argc, char *argv[]);

void main_usage(void);

int
main(int argc, char *argv[])
{
	char *mode = NULL;
	int val = 0;
	mode_t ret;

	if (argc < 2)
		main_usage();
	else
	{
		mode = strdup(argv[1]);
		ret = find_mode(mode);
		argc--;
		argv++;
		switch (ret)
		{
			case PARSE:
				val = parse_main(argc, argv);
				break;
			case TRIMEND:
				val = trimend_main(argc, argv);
				break;
			case PAIR:
				val = pair_main(argc, argv);
				break;
			case ERROR:
				val = 1;
				fprintf(stderr, "Mode \'%s\' not recognized.\n\n", mode);
				main_usage();
		}
	}

	/* Free allocated memory */
	free(mode);

	return val;
}

void
main_usage(void)
{
	fputs("Usage: ddradseq [MODE] [OPTIONS] [FASTQ] ...\n\n", stderr);
	fputs("Valid modes:\n", stderr);
	fputs("  parse	   Parses input fastQ by standard Illumina or custom adapter\n", stderr);
	fputs("  trimend  Trims the 3\' end of fastQ sequences\n", stderr);
	fputs("  pair	   Aligns mated pairs in two fastQ input files\n\n", stderr);
	fputs("Contact dgarriga@bioinformed.org for support.\n\n", stderr);
}
