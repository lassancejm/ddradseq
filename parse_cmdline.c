/* file: parse_cmdline.c
 * description: Function to read command line parameters
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>
#include <getopt.h>
#include <ctype.h>
#include <libgen.h>
#include "ddradseq.h"


CMD *
parse_cmdline(int argc, char *argv[], char *runmode)
{
	char *datec = NULL;
	int c = 0;
	runmode_t ret;
	size_t size = 0;
	time_t rawtime;
	struct tm * timeinfo;
	CMD *cp = NULL;

	/* Allocate memory for command line option structure */
	if ((cp = malloc(sizeof (CMD))) == NULL)
	{
		fputs("Error allocating memory for command line structure.\n", stderr);
		return NULL;
	}

	/* Initialize default values on command line data structure */
	cp->default_dir = false;
	cp->trim_barcode = false;
	cp->num_threads = 1;
	cp->parentdir = NULL;
	cp->outdir = NULL;
	cp->forfastq = NULL;
	cp->revfastq = NULL;
	cp->csvfile = NULL;

	/* Construct the date string for output directory */
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	datec = malloc(DATELEN + 2u);
	assert(datec != NULL);
	strftime(datec, DATELEN + 2u, "/%F/", timeinfo);

	while (1)
	{
		static struct option long_options[] =
		{
			{"trim",	no_argument,	   0, 't'},
			{"help",	no_argument,	   0, 'h'},
			{"out",		required_argument, 0, 'o'},
			{"csv",		required_argument, 0, 'c'},
			{"threads", required_argument, 0, 'n'},
			{0, 0, 0, 0}
		};

		int option_index = 0;

		c = getopt_long(argc, argv, "o:n:c:th", long_options, &option_index);
		if (c == -1) break;

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
				cp->parentdir = strdup(optarg);
				size = strlen(cp->parentdir);
				if (cp->parentdir[size - 1] == '/')
				{
					size += 12u;
					cp->outdir = malloc(size);
					strcpy(cp->outdir, cp->parentdir);
					strcat(cp->outdir, datec+1);
				}
				else
				{
					size += 13u;
					cp->outdir = malloc(size);
					strcpy(cp->outdir, cp->parentdir);
					strcat(cp->outdir, datec);
				}
				break;
			case 'c':
				cp->csvfile = strdup(optarg);
				break;
			case 'n':
				cp->num_threads = atoi(optarg);
				break;
			case 't':
				cp->trim_barcode = true;
				break;
			case 'h':
				free(cp);
				return NULL;
			case '?':
				free(cp);
				return NULL;
			default:
				free(cp);
				return NULL;
		}
	}

	/* Get run mode */
	ret = find_mode(runmode);

	/* Parse non-optioned arguments */
	if (ret == PARSE)
	{
		if ((optind + 2) > argc)
		{
			fputs("ERROR: two paired fastQ files are required as input\n\n",
				  stderr);
			free(cp);
			return NULL;
		}
		else if (cp->csvfile == NULL)
		{
			fputs("ERROR: \'--csv\' switch is mandatory in parse mode\n\n",
				  stderr);
			free(cp);
			return NULL;
		}
		else
		{
			cp->forfastq = strdup(argv[optind]);
			cp->revfastq = strdup(argv[optind + 1]);
			if (cp->parentdir == NULL)
			{
				char *fullpath = strdup(argv[optind]);
				cp->parentdir = dirname(fullpath);
				size = strlen(cp->parentdir);
				cp->outdir = malloc(size + 13u);
				strcpy(cp->outdir, cp->parentdir);
				strcat(cp->outdir, datec);
				cp->default_dir = true;
				free(fullpath);
			}
			return cp;
		}
	}
	else
	{
		if ((optind + 1) > argc)
		{
			fputs("ERROR: need the parent output directory as input\n\n",
				  stderr);
			free(cp);
			return NULL;
		}
		else
		{
			cp->outdir = strdup(argv[optind]);
			return cp;
		}
	}
}
