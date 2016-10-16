/* file: log_init.c
 * description: Initialize the ddradseq log file
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/utsname.h>
#include <sys/sysinfo.h>
#include "ddradseq.h"

int log_init(const CMD *cp)
{
	char *user = NULL;
	const double gigabyte = 1024 * 1024 * 1024;
	struct utsname host;
	struct sysinfo si;

	/* Print logfile header */
	fputs("***ddradseq LOG FILE***\n", lf);

	/* Get time information */
	get_timestr(&timestr[0]);

	/* Get user name */
	user = getenv("USER");

	/* Get host info */
	uname(&host);

	/* Get some system specs */
	sysinfo(&si);

	/* Print information on starting parameters */
	fprintf(lf, "[ddradseq: %s] INFO -- user specified directory %s "
	        "for input.\n", timestr, cp->parentdir);
	fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as database "
	        "file.\n", timestr, cp->csvfile);
	fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as output "
	        "directory.\n", timestr, cp->outdir);
	fprintf(lf, "[ddradseq: %s] INFO -- program will use edit distance of "
	        "%d base difference.\n", timestr, cp->dist);
	fprintf(lf, "[ddradseq: %s] INFO -- program has started in %s mode ",
	        timestr, cp->mode);
	if (user)
		fprintf(lf, "by user \'%s\' ", user);
	fprintf(lf, "on host \'%s\' (%s)", host.nodename, host.release);
	fputc('\n', lf);
	fprintf(lf, "[ddradseq: %s] INFO -- host has %5.1f Gb total RAM and %5.1f Gb "
	        "free RAM.\n", timestr, si.totalram/gigabyte, si.freeram/gigabyte);
	return 0;
}
