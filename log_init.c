/* file: log_init.c
 * description: Initialize the ddradseq log file
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/utsname.h>
#include <sys/sysinfo.h>
#include "ddradseq.h"

extern int errno;

int log_init(CMD *cp)
{
	char *logpath = NULL;
	const char *logfile = "/ddradseq.log";
	char *user = NULL;
	const double gigabyte = 1024 * 1024 * 1024;
	struct utsname host;
	struct sysinfo si;

	/* Get current working directory */
	logpath = get_current_dir_name();
	logpath = realloc(logpath, strlen(logpath) + strlen(logfile) + 1);
	strcat(logpath, logfile);

	/* Open log file for writing */
	cp->lf = fopen(logpath, "a");
	if (!cp->lf)
	{
		perror("Failed to open logfile");
		return 1;
	}
	free(logpath);

	/* Print logfile header */
	fputs("***ddradseq LOG FILE***\n", cp->lf);

	/* Get user name */
	user = getenv("USER");

	/* Get host info */
	uname(&host);

	/* Get some system specs */
	sysinfo(&si);

	/* Print information on starting parameters */
	loginfo(cp->lf, "user specified directory %s for input.\n", cp->parent_indir);
	loginfo(cp->lf, "user specified \'%s\' as database file.\n", cp->csvfile);
	loginfo(cp->lf, "user specified \'%s\' as output directory.\n", cp->parent_outdir);
	loginfo(cp->lf, "output will be written to \'%s\'.\n", cp->outdir);
	loginfo(cp->lf, "program will use edit distance of %d base difference.\n", cp->dist);
	loginfo(cp->lf, "program has started in %s mode ", cp->mode);
	if (user)
		loginfo(cp->lf, "by user \'%s\' ", user);
	loginfo(cp->lf, "on host \'%s\' (%s)\n", host.nodename, host.release);
	loginfo(cp->lf, "host has %5.1f Gb total RAM and %5.1f Gb free RAM.\n",
	        si.totalram/gigabyte, si.freeram/gigabyte);

	return 0;
}
