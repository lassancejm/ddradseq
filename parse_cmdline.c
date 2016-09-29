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


CMD *
parse_cmdline(int argc, char *argv[], char *runmode)
{
	char datestr[80];
	char *datec = NULL;
	int c = 0;
	runmode_t ret;
	size_t size = 0;
	time_t rawtime;
	struct tm * timeinfo;
	CMD *cp = NULL;
	FILE *lf;

	/* Open the log file */
	lf = fopen(logfile, "a");

	/* Allocate memory for command line option structure */
	if ((cp = malloc(sizeof (CMD))) == NULL)
	{
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(datestr, 80, "%c", timeinfo);
		fprintf(lf, "[ddradseq: %s] ERROR -- cannot allocate memory for command line structure.\n", datestr);
		fclose(lf);
		exit(EXIT_FAILURE);
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
		strftime(datestr, 80, "%c", timeinfo);
		fprintf(lf, "[ddradseq: %s] ERROR -- cannot allocate memory for datec.\n", datestr);
		fclose(lf);
		exit(EXIT_FAILURE);
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
				size = strlen(cp->parentdir);
				if (cp->parentdir[size - 1] == '/')
				{
					size += 12u;
					cp->outdir = malloc(size);
					strcpy(cp->outdir, cp->parentdir);
					strcat(cp->outdir, datec + 1);
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
			case 'd':
				cp->dist = atoi(optarg);
				break;
			case 'h':
				free(cp);
				free(datec);
				fclose(lf);
				return NULL;
			case '?':
				free(cp);
				free(datec);
				fclose(lf);
				return NULL;
			default:
				free(datec);
				free(cp);
				fclose(lf);
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
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(datestr, 80, "%c", timeinfo);
			fprintf(lf, "[ddradseq: %s] ERROR -- two paired fastQ files are required as input.\n", datestr);
			free(cp);
			free(datec);
			fclose(lf);
			return NULL;
		}
		else if (cp->csvfile == NULL)
		{
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(datestr, 80, "%c", timeinfo);
			fprintf(lf, "[ddradseq: %s] ERROR -- \'--csv\' switch is mandatory in parse mode.\n", datestr);
			free(cp);
			free(datec);
			fclose(lf);
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
				cp->default_dir = 1;
				free(fullpath);
			}
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(datestr, 80, "%c", timeinfo);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as input fastQ.1.\n", datestr, cp->forfastq);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as input fastQ.2.\n", datestr, cp->revfastq);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as database file.\n", datestr, cp->csvfile);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified \'%s\' as output directory.\n", datestr, cp->outdir);
			fprintf(lf, "[ddradseq: %s] INFO -- program will use edit distance =  %d.\n", datestr, cp->dist);
			free(datec);
			fclose(lf);
			return cp;
		}
	}
	else
	{
		if ((optind + 1) > argc)
		{
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(datestr, 80, "%c", timeinfo);
			fprintf(lf, "[ddradseq: %s] ERROR -- need the parent output directory as input.\n", datestr);
			free(datec);
			free(cp);
			fclose(lf);
			return NULL;
		}
		else
		{
			cp->outdir = strdup(argv[optind]);
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			strftime(datestr, 80, "%c", timeinfo);
			fprintf(lf, "[ddradseq: %s] INFO -- user specified directory %s for input.\n", datestr, cp->outdir);
			free(datec);
			fclose(lf);
			return cp;
		}
	}
}
