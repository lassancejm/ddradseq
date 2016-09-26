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
#include <unistd.h>
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

	opterr = 0;

	/* Allocate memory for command line option structure */
	if ((cp = malloc(sizeof (CMD))) == NULL)
	{
		fputs("Error allocating memory for command line structure.\n", stderr);
		return NULL;
	}


	/* Initialize default values on command line data structure */
	cp->default_dir = false;
	cp->is_reverse = false;
	cp->num_threads = 1;
	cp->parentdir = NULL;
	cp->outdir = NULL;
	cp->filename = NULL;
	cp->csvfile = NULL;

	while ((c = getopt (argc, argv, "o:t:i:rh")) != -1)
	{
		switch (c)
		{
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
			case 'r':
				cp->is_reverse = true;
				break;
			case 'i':
				cp->csvfile = strdup(optarg);
				break;
			case 't':
				cp->num_threads = atoi(optarg);
				break;
			case 'h':
				free(cp);
				return NULL;
			case '?':
				if (optopt == 'c')
					fprintf(stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint(optopt))
					fprintf(stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
				free(cp);
				return NULL;
			default:
				free(cp);
				return NULL;
		}
	}

	if ((optind + 1) > argc)
	{
		free(cp);
		return NULL;
	}
	else
	{
		cp->filename = strdup(argv[optind]);;
		if (cp->parentdir == NULL)
		{
			char *fullpath = strdup(argv[optind]);
			cp->parentdir = dirname(fullpath);
			size = strlen(cp->parentdir);
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
			cp->default_dir = true;
		}
		return cp;
	}
}
