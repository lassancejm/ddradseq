/* file: ddradseq.c
 * description: Entry point for the ddradseq program
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "ddradseq.h"

/* Function prototypes */
extern int parse_main(const CMD*);
extern int trimend_main(const CMD*);
extern int pair_main(const CMD*);

extern int errno;

int main(int argc, char *argv[])
{
	char *cwd = NULL;
	char logfname[] = "/ddradseq.log\0";
	int ret = 0;
	CMD *cp = NULL;

	/* Parse the command line options */
	cp = parse_cmdline(argc, argv);
	if (!cp)
		return 1;

	/* Get current working directory */
	cwd = getcwd(logfile, sizeof(logfile));
	if (!cwd)
	{
		perror("Failed to get current working directory");
		return 1;
	}
	strcat(logfile, logfname);

	/* Open log file for writing */
	lf = fopen(logfile, "a");
	if (!lf)
	{
		perror("Failed to open logfile");
		return 1;
	}

	/* Initialize the log files */
	ret = log_init(cp);
	if (ret)
		return 1;

	/* Run the parse pipeline stage */
	if (string_equal(cp->mode, "parse") || string_equal(cp->mode, "all"))
	{
		ret = parse_main(cp);
		if (ret)
			return 1;
	}

	/* Run the pair pipeline stage */
	if (string_equal(cp->mode, "pair") || string_equal(cp->mode, "all"))
	{
		ret = pair_main(cp);
		if (ret)
			return 1;
	}

	/* Run the trimend pipeline stage */
	if (string_equal(cp->mode, "trimend") || string_equal(cp->mode, "all"))
	{
		ret = trimend_main(cp);
		if (ret)
			return 1;
	}

	/* Close logfile output stream */
	fclose(lf);

	/* Free memory for command line data structure from heap */
	free_cmdline(cp);

	return 0;
}
