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
extern int parse_main(int argc, char *argv[]);
extern int trimend_main(int argc, char *argv[]);
extern int pair_main(int argc, char *argv[]);
void main_usage(void);

/* Function pointer */
typedef int (*func_t)( int x, char *y[]);

int
main(int argc, char *argv[])
{
	char logfname[] = "/ddradseq.log\0";
	char *runmode = NULL;
	char version_str[] = "ddradseq v0.9-alpha";
	int val = EXIT_SUCCESS;
	runmode_t ret;
	func_t f[3];

	f[0] = parse_main;
	f[1] = trimend_main;
	f[2] = pair_main;

	/* Get current working directory */
	if (getcwd(logfile, sizeof(logfile)) == NULL)
	{
		fputs("ERROR: Failed to get current working directory.\n", stderr);
		return EXIT_FAILURE;
	}
	strcat(logfile, logfname);

	/* Open log file for writing */
	lf = fopen(logfile, "a");
	
	if (argc < 2)
		main_usage();
	else
	{
		if (strcmp(argv[1], "--version") == 0 ||
			strcmp(argv[1], "-v") == 0)
		{
			fprintf(stdout, "%s\n", version_str);
			return EXIT_SUCCESS;
		}
		runmode = strdup(argv[1]);
		ret = find_mode(runmode);
		if (ret == ERROR)
		{
			fprintf(stderr, "ERROR: Mode \'%s\' not recognized.\n\n", runmode);
			main_usage();
			val = EXIT_FAILURE;
		}
		else
		{
			log_init(ret);
			argc--;
			argv++;
			val = f[ret](argc, argv);
		}
	}

	/* Close logfile output stream */
	fclose(lf);

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
