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
#include "ddradseq.h"

/* Function prototypes */
extern int parse_main(int argc, char *argv[]);
extern int trimend_main(int argc, char *argv[]);
extern int pair_main(int argc, char *argv[]);
void main_usage(void);

/* Function pointer */
typedef int (*func_t)( int x, char *y[]);

int
main(int argc, char *argv[])
{
	char *runmode = NULL;
	char version_str[] = "ddradseq v0.9-alpha";
	int val = 0;
	runmode_t ret;
	func_t f[3];

	f[0] = parse_main;
	f[1] = trimend_main;
	f[2] = pair_main;

	if (argc < 2)
		main_usage();
	else
	{
		if (strcmp(argv[1], "--version") == 0 ||
			strcmp(argv[1], "-v") == 0)
		{
			fprintf(stdout, "%s\n", version_str);
			return 0;
		}
		runmode = strdup(argv[1]);
		ret = find_mode(runmode);
		if (ret == ERROR)
		{
			val = 1;
			fprintf(stderr, "Mode \'%s\' not recognized.\n\n", runmode);
			main_usage();
		}
		else
		{
			argc--;
			argv++;
			val = f[ret](argc, argv);
		}
	}

	/* Free allocated memory */
	free(runmode);

	return val;
}

void
main_usage(void)
{
	fputs("Usage: ddradseq [MODE] [OPTIONS] [INPUT FILES/DIRECTORY]\n\n", stderr);
	fputs("Valid modes are:\n", stderr);
	fputs("	 parse	   Parses input fastQ by standard Illumina and custom adapter\n", stderr);
	fputs("	 pair	   Aligns mated pairs in two fastQ input files\n", stderr);
	fputs("	 trimend   Trims the 3\' end of fastQ reverse sequences\n\n", stderr);
	fputs("The modes are intended to be run in the above order.\n", stderr);
	fputs("Use \'ddradseq -v\' or \'ddradseq --version\' to see software version\n", stderr);
	fputs("For development information, see https://github.com/dgarriga/ddradseq\n", stderr);
	fputs("Contact dgarriga@lummei.net for support.\n\n", stderr);
}
