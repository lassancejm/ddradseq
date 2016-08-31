/* file: parse_pool.c
 * description: Worker program for the ddRadSeq pipeline to sort fastQ
 * entries by Illumina index and write each pool to separate file
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: August 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <libgen.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fasta.h>

/* Define some main return errors */
#define INCORRECT_PARAMETERS 2
#define DIR_NOT_WRITABLE 3
#define CREATE_DIR_FAILED 4
#define DATABASE_CREATION_FAILED 5
#define ADAPTER_DB_FAILED 6
#define PARSE_FAILED 7
#define ALLOC_ERROR 8

/* Define data structures */
struct cmdparam
{
    int default_dir;
	int num_threads;
    char *parentdir;
    char *outdir;
    char *filename;
    char *csvfile;
};

/* Function prototypes */
struct cmdparam *parse_cmdline (int argc, char *argv[]);
void usage (char *program);

int
main (int argc, char *argv[])
{
    int writable = 0;
    int status = 0;
    struct cmdparam *cp = NULL;
    FQDB *reads = NULL;
    STDDB *adapt = NULL;

    opterr = 0;

    /* Parse the command line options */
    if ((cp = parse_cmdline (argc, argv)) == NULL)
        {
            usage (argv[0]);
            return INCORRECT_PARAMETERS;
        }

    /* Check if parent output directory is writable */
    writable = access (cp->parentdir, W_OK);
    if (writable != 0)
        {
            fprintf (stderr, "Error: cannot write to directory %s.\n",
                    cp->parentdir);
            return DIR_NOT_WRITABLE;
        }
    else
        {
			/* Test if output directory already exists */
			DIR *d = opendir (cp->outdir);

			/* If the directory doesn't already exist, create it */
			if (d)
				{
					closedir (d);
				}
			else if (ENOENT == errno)
				{
					status = mkdir (cp->outdir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
					if (status < 0)
						{
							fprintf (stderr, "Error: failed to create output "
									 "directory %s.\n", cp->outdir);
							return CREATE_DIR_FAILED;
						}
				}
        }

    /* Read the fastQ input files */
    reads = read_fastq (1, cp->filename);
    if (reads == NULL)
		{
			fprintf (stderr, "Error: cannot read input file %s.\n",
					 cp->filename);
			return DATABASE_CREATION_FAILED;
		}

    /* Read the adapter CSV input file and create database */
    adapt = create_adapter_database (cp->csvfile);
    if (adapt == NULL)
		{
			fprintf (stderr, "Error cannot read CSV file %s.\n", cp->csvfile);
			return ADAPTER_DB_FAILED;
		}

    /* Parse the reads by Illumina index and write to separate files */
    if ((parse_fastq_illumina (reads, adapt, cp->outdir, cp->num_threads)) < 0)
        {
            fputs ("Error in parse_fastq_illumina.\n", stderr);
            return PARSE_FAILED;
        }

    /* Deallocate memory */
	if (!cp->default_dir)
		{
			free (cp->parentdir);
		}
	free (cp->outdir);
    free (cp->filename);
    free (cp->csvfile);
    free (cp);

    return EXIT_SUCCESS;
}

struct cmdparam *
parse_cmdline (int argc, char *argv[])
{
    int c = 0;
    size_t size = 0;
    struct cmdparam *cp = NULL;

    /* Allocate memory for command line option structure */
    if ((cp = malloc (sizeof (struct cmdparam))) == NULL)
        {
            fprintf (stderr, "[parse_pool:%s:%d] Error allocating memory "
                     "for command line structure.\n", __func__, __LINE__);
            exit (ALLOC_ERROR);
        }

	/* Initialize default values on command line data structure */
	cp->default_dir = 0;
	cp->num_threads = 1;
    cp->parentdir = NULL;
    cp->outdir = NULL;
    cp->filename = NULL;
    cp->csvfile = NULL;

    while ((c = getopt (argc, argv, "o:t:h")) != -1)
        {
            switch (c)
                {
                    case 'o':
                        cp->parentdir = strdup (optarg);
                        size = strlen (cp->parentdir);
                        if (cp->parentdir[size - 1] == '/')
                            {
                                cp->outdir = malloc (size + 6);
                                strcpy (cp->outdir, cp->parentdir);
                                strcat (cp->outdir, "pool/");
                            }
                        else
                            {
                                cp->outdir = malloc (size + 7);
                                strcpy (cp->outdir, cp->parentdir);
                                strcat (cp->outdir, "/pool/");
                            }
						break;
					case 't':
						cp->num_threads = atoi (optarg);
                        break;
                    case 'h':
                        free (cp);
                        return NULL;
                    case '?':
                        if (optopt == 'c')
                            {
                                fprintf (stderr, "Option -%c requires "
                                         "an argument.\n", optopt);
                            }
                        else if (isprint (optopt))
                            {
                                fprintf (stderr, "Unknown option `-%c'."
                                         "\n", optopt);
                            }
                        else
                            {
                                fprintf (stderr, "Unknown option "
                                         "character `\\x%x'.\n", optopt);
                            }
                        free (cp);
                        return NULL;
                    default:
                        free (cp);
                        return NULL;
                }
        }

    if ((optind + 2) > argc)
        {
            free (cp);
            return NULL;
        }
    else
        {
            cp->filename = strdup (argv[optind]);;
            cp->csvfile = strdup (argv[optind + 1]);
            if (cp->parentdir == NULL)
                {
                    char *fullpath = strdup (argv[optind]);
                    cp->parentdir = dirname (fullpath);
                    size = strlen (cp->parentdir);
					cp->outdir = malloc (size + 7);
					strcpy (cp->outdir, cp->parentdir);
					strcat (cp->outdir, "/pool/");
					cp->default_dir = 1;
                }
            return cp;
        }
}

void
usage (char *program)
{
    fprintf (stderr, "Usage : %s [OPTIONS] [FASTQ] [CSV]\n", program);
    fputs ("Parse fastQ file into separate files by Illumina index\n\n",
           stderr);
    fputs ("Available options\n", stderr);
    fputs ("  -o  DIR    Parent directory to write output files "
           "[default: same as input fastQ]\n", stderr);
	fputs ("  -t  INT    Number of threads for concurrency [default: 1]\n", stderr);
    fputs ("  -h         Display this help message\n\n", stderr);
}

