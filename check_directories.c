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
	int writable = 0;
	int status = 0;
	size_t strl = 0;
	khint_t i = 0;
	khint_t j = 0;
	khash_t(pool) *p = NULL;
	POOL *pl = NULL;
	DIR *d;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Check if parent output directory is writable */
	writable = access(cp->parentdir, W_OK);
	if (writable != 0)
	{
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Cannot write to directory "
		        "\'%s\'.\n", timestr, __func__, __LINE__, cp->parentdir);
		return 1;
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
				fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to create output "
				        "directory \'%s\'.\n", timestr, __func__, __LINE__, cp->outdir);
				return 1;
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
						fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to create "
						        "flowcell-level output directory \'%s\'.\n", timestr,
						        __func__, __LINE__, flowdir);
						return 1;
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
								fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to create "
								        "pool-level directory \'%s\'.\n", timestr, __func__,
								        __LINE__, pooldir);
								return 1;
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
								fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to create "
								        "output directory \'%s\'.\n", timestr, __func__,
								        __LINE__, parsedir);
								return 1;
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
								fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to create "
								        "output directory \'%s\'.\n", timestr, __func__,
								        __LINE__, pairdir);
								return 1;
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
								fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Failed to create "
								        "output directory \'%s\'.\n", timestr, __func__,
								        __LINE__, trimdir);
								return 1;
							}
						}
						free(trimdir);
					}
				}
			}
		}
	}

	return 0;
}
