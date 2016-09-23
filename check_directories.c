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
	int writable = 0;
	int status = 0;
	size_t strl = 0;
	khint_t i = 0;
	khint_t j = 0;
	khash_t(pool) *p = NULL;
	POOL *pl = NULL;
	DIR *d;

	/* Check if parent output directory is writable */
	writable = access(cp->parentdir, W_OK);
	if (writable != 0)
	{
		fprintf(stderr, "Error: cannot write to directory "
		        "\'%s\'.\n", cp->parentdir);
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
				fprintf(stderr, "Error: failed to create output "
				        "directory \'%s\'.\n", cp->outdir);
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
					status = mkdir(flowdir,  S_IRWXU | S_IRGRP | S_IXGRP |
			                                 S_IROTH | S_IXOTH);
					if (status < 0)
					{
						fprintf(stderr, "Error: failed to create output "
						        "directory \'%s\'.\n", flowdir);
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
						pooldir = pl->dirpath;
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
								fprintf(stderr, "Error: failed to create "
								        "output directory \'%s\'.\n", pooldir);
								return 1;
				 			}
				 		}
					}
				}
			}
		}
 	}
	return 0;
}
