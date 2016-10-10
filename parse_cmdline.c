/* file: parse_cmdline.c
 * description: Populates command line data structure from program arguments
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <libgen.h>
#include "ddradseq.h"

CMD *parse_cmdline(int argc, char *argv[])
{
	char *datec = NULL;
	char *tmpdir = NULL;
	char version_str[] = "ddradseq v1.0-beta";
	int c = 0;
	int r = 0;
	size_t strl = 0;
	time_t rawtime;
	struct tm *timeinfo;
	CMD *cp = NULL;

	/* Allocate memory for command line option structure */
	cp = malloc(sizeof(CMD));
	if (UNLIKELY(!cp))
	{
		error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return NULL;
	}

	/* Initialize default values on command line data structure */
	cp->parentdir = NULL;
	cp->outdir = NULL;
	cp->csvfile = NULL;
	cp->mode = NULL;
	cp->dist = 1;
	cp->score = 100;
	cp->gapo = 5;
	cp->gape = 2;

	while (1)
	{
		static struct option long_options[] =
		{
			{"help",    no_argument,       0, 'h'},
			{"version", no_argument,       0, 'v'},
			{"mode",    required_argument, 0, 'm'},
			{"dist",    required_argument, 0, 'd'},
			{"out",	    required_argument, 0, 'o'},
			{"score",   required_argument, 0, 's'},
			{"gapo",    required_argument, 0, 'g'},
			{"gape",    required_argument, 0, 'e'},
			{"csv",	    required_argument, 0, 'c'},
			{0, 0, 0, 0}
		};

		int option_index = 0;

		c = getopt_long(argc, argv, "o:d:c:m:s:hv", long_options, &option_index);
		if (c == -1)
			break;

		switch (c)
		{
			case 0:
				if (long_options[option_index].flag != 0)
					break;
				printf("option %s", long_options[option_index].name);
				if (optarg)
					printf(" with arg %s", optarg);
				printf("\n");
				break;
			case 'o':
				/* Construct the date string for output directory */
				datec = malloc(DATELEN + 3u);
				if (UNLIKELY(!datec))
				{
					error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return NULL;
				}
				time(&rawtime);
				timeinfo = localtime(&rawtime);
				strftime(datec, DATELEN + 3u, "/ddradseq-%F/", timeinfo);
				tmpdir = strdup(optarg);
				strl = strlen(tmpdir);
				if (tmpdir[strl - 1] == '/')
				{
					strl += 22u;
					cp->outdir = malloc(strl);
					if (UNLIKELY(!cp->outdir))
					{
						error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return NULL;
					}
					strcpy(cp->outdir, tmpdir);
					strcat(cp->outdir, datec + 1);
				}
				else
				{
					strl += 23u;
					cp->outdir = malloc(strl);
					if (UNLIKELY(!cp->outdir))
					{
						error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return NULL;
					}
					strcpy(cp->outdir, tmpdir);
					strcat(cp->outdir, datec);
				}
				free(tmpdir);
				free(datec);
				break;
			case 'm':
				cp->mode = strdup(optarg);
				if ((r = strcmp(cp->mode, "parse")) != 0 ||
				    (r = strcmp(cp->mode, "pair")) != 0 ||
				    (r = strcmp(cp->mode, "trimend")) != 0)
				{
					fprintf(stderr, "ERROR: %s is not a valid mode.\n", cp->mode);
					return NULL;
				}
				break;
			case 's':
				cp->score = atoi(optarg);
				break;
			case 'g':
				cp->gapo = atoi(optarg);
				break;
			case 'e':
				cp->gape = atoi(optarg);
				break;
			case 'c':
				cp->csvfile = strdup(optarg);
				break;
			case 'd':
				cp->dist = atoi(optarg);
				break;
			case 'v':
				fprintf(stdout, "%s\n", version_str);
				exit(EXIT_SUCCESS);
			case 'h':
				return NULL;
			case '?':
				return NULL;
			default:
				return NULL;
		}
	}

	/* Do sanity check on command line options */
	if (!cp->mode)
	{
		cp->mode = malloc(4u);
		strcpy(cp->mode, "all");
	}
	if (!cp->csvfile && ((r = strcmp(cp->mode, "parse")) == 0 || (r = strcmp(cp->mode, "all")) == 0))
	{
		fputs("ERROR: \'--csv\' switch is mandatory when running parse mode.\n", stderr);
		return NULL;
	}
	if (!cp->outdir && ((r = strcmp(cp->mode, "parse")) == 0 || (r = strcmp(cp->mode, "all")) == 0))
	{
		fputs("ERROR: \'--out\' switch is mandatory when running parse mode.\n", stderr);
		return NULL;
	}

	/* Parse non-optioned arguments */
	if (optind + 1 > argc)
	{
		fputs("ERROR: need the fastQ parent directory as input.\n", stderr);
		return NULL;
	}
	else
	{
		cp->parentdir = strdup(argv[optind]);
		return cp;
	}
}
