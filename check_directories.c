/* file: check_directories.c
 * description: Function to check write permissions on parent directory
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <dirent.h>
#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "khash.h"
#include "ddradseq.h"

int
check_directories(CMD *cp, khash_t(pool_hash) *h)
{
	char *pooldir = NULL;
	char *flowdir = NULL;
	char *parsedir = NULL;
	char *pairdir = NULL;
	char *trimdir = NULL;
	char datestr[80];
	int writable = 0;
	int status = 0;
	size_t strl = 0;
	time_t rawtime;
	struct tm * timeinfo;
	khint_t i = 0;
	khint_t j = 0;
	khash_t(pool) *p = NULL;
	POOL *pl = NULL;
	DIR *d;
	FILE *lf;

	/* Open logfile for writing */
	lf = fopen(logfile, "a");

	/* Get time info */
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(datestr, 80, "%c", timeinfo);

	/* Check if parent output directory is writable */
	writable = access(cp->parentdir, W_OK);
	if (writable != 0)
	{
		fprintf(lf, "[ddradseq: %s] ERROR -- cannot write to directory \'%s\'.\n", datestr, cp->parentdir);
		fclose(lf);
		exit(EXIT_FAILURE);
	}
	else
	{
		/* Test if output directory already exists */
		d = opendir(cp->outdir);

		/* If the directory doesn't already exist, create it */
		if (d)
			closedir(d);
		else if (ENOENT == errno)
		{
			status = mkdir(cp->outdir, S_IRWXU | S_IRGRP | S_IXGRP |
									   S_IROTH | S_IXOTH);
			if (status < 0)
			{
				fprintf(lf, "[ddradseq: %s] ERROR -- failed to create output directory \'%s\'.\n", datestr, cp->outdir);
				fclose(lf);
				exit(EXIT_FAILURE);
			}
		}

		/* Check and create pool subdirectories */
		for (i = kh_begin(h); i != kh_end(h); i++)
		{
			if (kh_exist(h, i))
			{
				strl = strlen(cp->outdir) + strlen(kh_key(h, i));
				flowdir = malloc(strl + 1u);
				strcpy(flowdir, cp->outdir);
				strcat(flowdir, kh_key(h, i));
				d = opendir(flowdir);
				if (d)
					closedir(d);
				else if (ENOENT == errno)
				{
					status = mkdir(flowdir, S_IRWXU | S_IRGRP | S_IXGRP |
											S_IROTH | S_IXOTH);
					if (status < 0)
					{
						fprintf(lf, "[ddradseq: %s] ERROR -- failed to create flowcell-level output directory \'%s\'.\n", datestr, flowdir);
						fclose(lf);
						exit(EXIT_FAILURE);
					}
				}
				free(flowdir);
				p = kh_value(h, i);
				for (j = kh_begin(p); j != kh_end(p); j++)
				{
					if (kh_exist(p, j))
					{
						pl = kh_value(p, j);
						pooldir = pl->poolpath;
						d = opendir(pooldir);

						/* If the subsubdirectory doesn't already exist */
						/* create it */
						if (d)
							closedir(d);
						else if (ENOENT == errno)
						{
							status = mkdir(pooldir, S_IRWXU | S_IRGRP |
													S_IXGRP | S_IROTH | S_IXOTH);
							if (status < 0)
							{
								fprintf(lf, "[ddradseq: %s] ERROR -- failed to create pool-level directory \'%s\'.\n", datestr, pooldir);
								fclose(lf);
								exit(EXIT_FAILURE);
							}
						}
						strl = strlen(pooldir);
						parsedir = malloc(strl + 7u);
						strcpy(parsedir, pooldir);
						strcat(parsedir, "/parse");
						d = opendir(parsedir);
						/* If the subsubdirectory doesn't already exist */
						/* create it */
						if (d)
							closedir(d);
						else if (ENOENT == errno)
						{
							status = mkdir(parsedir, S_IRWXU | S_IRGRP |
													 S_IXGRP | S_IROTH | S_IXOTH);
							if (status < 0)
							{
								fprintf(lf, "[ddradseq: %s] ERROR -- failed to create output directory \'%s\'.\n", datestr, parsedir);
								fclose(lf);
								exit(EXIT_FAILURE);
							}
						}
						free(parsedir);
						strl = strlen(pooldir);
						pairdir = malloc(strl + 7u);
						strcpy(pairdir, pooldir);
						strcat(pairdir, "/pairs");
						d = opendir(pairdir);
						/* If the subsubdirectory doesn't already exist */
						/* create it */
						if (d)
							closedir(d);
						else if (ENOENT == errno)
						{
							status = mkdir(pairdir, S_IRWXU | S_IRGRP |
													S_IXGRP | S_IROTH | S_IXOTH);
							if (status < 0)
							{
								fprintf(lf, "[ddradseq: %s] ERROR -- failed to create output directory \'%s\'.\n", datestr, pairdir);
								fclose(lf);
								exit(EXIT_FAILURE);
							}
						}
						free(pairdir);
						strl = strlen(pooldir);
						trimdir = malloc(strl + 7u);
						strcpy(trimdir, pooldir);
						strcat(trimdir, "/trime");
						d = opendir(trimdir);
						/* If the subsubdirectory doesn't already exist */
						/* create it */
						if (d)
							closedir(d);
						else if (ENOENT == errno)
						{
							status = mkdir(trimdir, S_IRWXU | S_IRGRP |
													S_IXGRP | S_IROTH | S_IXOTH);
							if (status < 0)
							{
								fprintf(lf, "[ddradseq: %s] ERROR -- failed to create output directory \'%s\'.\n", datestr, trimdir);
								fclose(lf);
								exit(EXIT_FAILURE);
							}
						}
						free(trimdir);
					}
				}
			}
		}
	}

	/* Close log file */
	fclose(lf);

	return 0;
}
