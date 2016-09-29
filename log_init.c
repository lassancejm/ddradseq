/* file: log_init.c
 * description: Initialize the ddradseq log file
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/utsname.h>
#include <sys/sysinfo.h>
#include "ddradseq.h"

void log_init(runmode_t runmode)
{
	char *user = NULL;
	char datestr[80];
	const double gigabyte = 1024 * 1024 * 1024;
	time_t rawtime;
	struct tm * timeinfo;
	struct utsname host;
	struct sysinfo si;
	FILE *lf;

	lf = fopen(logfile, "a");

	fputs("***ddradseq LOG FILE***\n", lf);

	/* Get time information */
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(datestr, 80, "%c", timeinfo);

	/* Get user name */
	user = getenv("USER");

	/* Get host info */
    uname(&host);

	/* Get some system specs */
	sysinfo (&si);

	switch(runmode)
	{
		case PARSE:
			fprintf(lf, "[ddradseq: %s] INFO -- program has started in parse mode ", datestr);
			break;
		case PAIR:
			fprintf(lf, "[ddradseq: %s] INFO -- program has started in pair mode ", datestr);
			break;
		case TRIMEND:
			fprintf(lf, "[ddradseq: %s] INFO -- program has started in trimend mode ", datestr);
			break;
		case ERROR:
			fprintf(lf, "[ddradseq: %s] ERROR -- program has started in unknown mode ", datestr);
			break;
	}
	if (user) fprintf(lf, "by user %s ", user);
	fprintf(lf, "on host %s %s", host.nodename, host.release);
	fputc('\n', lf);
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(datestr, 80, "%c", timeinfo);
	fprintf(lf, "[ddradseq: %s] INFO -- host has %5.1f Gb total RAM and %5.1f Gb free RAM.\n", datestr, si.totalram/gigabyte, si.freeram/gigabyte);
	fclose(lf);
}