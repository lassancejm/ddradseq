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
#include <time.h>
#include <getopt.h>
#include <libgen.h>
#include "ddradseq.h"

void parse_cmdline_deallocate(CMD *, char *);

CMD *
parse_cmdline(int argc, char *argv[], char *runmode)
{
	char *datec = NULL;
	int c = 0;
	runmode_t ret;
	size_t strl = 0;
	time_t rawtime;
	struct tm *timeinfo;
	CMD *cp = NULL;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Allocate memory for command line option structure */
	if ((cp = malloc(sizeof(CMD))) == NULL)
	{
		fputs("ERROR: Memory allocation failure.\n", stderr);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
		        timestr, __func__, __LINE__);
		return NULL;
	}

	/* Initialize default values on command line data structure */
	cp->default_dir = 0;
	cp->parentdir = NULL;
	cp->outdir = NULL;
	cp->forfastq = NULL;
	cp->revfastq = NULL;
	cp->csvfile = NULL;
	cp->dist = 1;

	/* Construct the date string for output directory */
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	if ((datec = malloc(DATELEN + 2u)) == NULL)
	{
		fputs("ERROR: Memory allocation failure.\n", stderr);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
		        timestr, __func__, __LINE__);
		return NULL;
	}
	strftime(datec, DATELEN + 2u, "/%F/", timeinfo);

	while (1)
	{
		static struct option long_options[] =
		{
			{"help",	no_argument,	   0, 'h'},
			{"dist",    required_argument, 0, 'd'},
			{"out",		required_argument, 0, 'o'},
			{"csv",		required_argument, 0, 'c'},
			{0, 0, 0, 0}
		};

		int option_index = 0;

		c = getopt_long(argc, argv, "o:d:c:h", long_options, &option_index);
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
				strl = strlen(cp->parentdir);
				if (cp->parentdir[strl - 1] == '/')
				{
					strl += 12u;
					if ((cp->outdir = malloc(strl)) == NULL)
					{
						fputs("ERROR: Memory allocation failure.\n", stderr);
						fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
						        timestr, __func__, __LINE__);
						return NULL;
					}
					strcpy(cp->outdir, cp->parentdir);
					strcat(cp->outdir, datec + 1);
				}
				else
				{
					strl += 13u;
					if ((cp->outdir = malloc(strl)) == NULL)
					{
						fputs("ERROR: Memory allocation failure.\n", stderr);
						fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
						        timestr, __func__, __LINE__);
						return NULL;
					}
					strcpy(cp->outdir, cp->parentdir);
					strcat(cp->outdir, datec);
				}
				break;
			case 'c':
				cp->csvfile = strdup(optarg);
				break;
			case 'd':
				cp->dist = atoi(optarg);
				break;
			case 'h':
				parse_cmdline_deallocate(cp, datec);
				return NULL;
			case '?':
				parse_cmdline_deallocate(cp, datec);
				return NULL;
			default:
				parse_cmdline_deallocate(cp, datec);
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
			fputs("ERROR: Paired fastQ files are required as input.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Paired fastQ files are "
			        "required as input.\n", timestr, __func__, __LINE__);
			parse_cmdline_deallocate(cp, datec);
			return NULL;
		}
		else if (cp->csvfile == NULL)
		{
			fputs("ERROR: \'--csv\' switch is mandatory in parse mode.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d \'--csv\' switch is mandatory in "
			        "parse mode.\n", timestr, __func__, __LINE__);
			parse_cmdline_deallocate(cp, datec);
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
				strl = strlen(cp->parentdir);
				cp->outdir = malloc(strl + 13u);
				strcpy(cp->outdir, cp->parentdir);
				strcat(cp->outdir, datec);
				cp->default_dir = 1;
				free(fullpath);
			}
			fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as input "
			        "fastQ.1.\n", timestr, cp->forfastq);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as input "
			        "fastQ.2.\n", timestr, cp->revfastq);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as database "
			        "file.\n", timestr, cp->csvfile);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as output "
			        "directory.\n", timestr, cp->outdir);
			fprintf(lf, "[ddradseq: %s] INFO -- program will use edit distance of "
			        "%d base difference.\n", timestr, cp->dist);
			parse_cmdline_deallocate(NULL, datec);
			return cp;
		}
	}
	else
	{
		if ((optind + 1) > argc)
		{
			fputs("ERROR: need the parent output directory as input.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Need the parent output "
			        "directory as input.\n", timestr, __func__, __LINE__);
			parse_cmdline_deallocate(cp, datec);
			return NULL;
		}
		else
		{
			cp->outdir = strdup(argv[optind]);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified directory %s "
			        "for input.\n", timestr, cp->outdir);
			parse_cmdline_deallocate(NULL, datec);
			return cp;
		}
	}
}

void
parse_cmdline_deallocate(CMD *cp, char *s)
{
	if (s) free(s);
	if (cp) free_cmdline(cp);	
}
