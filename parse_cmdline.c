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
#include <getopt.h>
#include <ctype.h>
#include <libgen.h>
#include "ddradseq.h"


CMD *
parse_cmdline(int argc, char *argv[], char *runmode)
{
	int c = 0;
	runmode_t ret;
	size_t size = 0;
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

	while (1)
	{
		static struct option long_options[] =
		{
			{"trim",    no_argument,       0, 't'},
			{"help",    no_argument,       0, 'h'},
			{"out",     required_argument, 0, 'o'},
			{"csv",     required_argument, 0, 'c'},
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
					cp->outdir = malloc(size + 9u);
					strcpy(cp->outdir, cp->parentdir);
                    ret = find_mode(runmode);
					switch (ret)
                    {
                        case PARSE:
                            strcat(cp->outdir, "parse/");
                            break;
                        case TRIMEND:
                            strcat(cp->outdir, "trimend/");
                            break;
                        case PAIR:
                            strcat(cp->outdir, "pair/");
                            break;
                        case ERROR:
                            return NULL;
                    }
				}
				else
				{
					cp->outdir = malloc(size + 10u);
					strcpy(cp->outdir, cp->parentdir);
                    ret = find_mode(runmode);
					switch (ret)
                    {
                        case PARSE:
                            strcat(cp->outdir, "/parse/");
                            break;
                        case TRIMEND:
                            strcat(cp->outdir, "/trimend/");
                            break;
                        case PAIR:
                            strcat(cp->outdir, "/pair/");
                            break;
                        case ERROR:
                            return NULL;
                    }
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

	if ((optind + 2) > argc)
	{
		fputs("ERROR: two paired fastQ files are required as input\n\n", stderr);
		free(cp);
		return NULL;
	}
	else if ((ret = find_mode(runmode)) == PARSE && cp->csvfile == NULL)
	{
		fputs("ERROR: \'--csv\' switch is mandatory in parse mode\n\n", stderr);
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
			cp->outdir = malloc(size + 10u);
			strcpy(cp->outdir, cp->parentdir);
            /*ret = find_mode(runmode);*/
            switch (ret)
            {
                case PARSE:
                    strcat(cp->outdir, "/parse/");
                    break;
                case TRIMEND:
                    strcat(cp->outdir, "/trimend/");
                    break;
                case PAIR:
                    strcat(cp->outdir, "/pair/");
                    break;
                case ERROR:
                    return NULL;
            }
			cp->default_dir = true;
		}
		return cp;
	}
}
