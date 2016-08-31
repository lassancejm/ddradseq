/* file: trim_3prime.c
 * description: Worker program for the ddRadSeq pipeline to trim adapter
 * sequence from 3' end of reverse reads
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: August 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fasta.h>

/* Define some main return errors */
#define INCORRECT_PARAMETERS 2
#define DIR_NOT_WRITABLE 3
#define CREATE_DIR_FAILED 4
#define DATABASE_CREATION_FAILED 5
#define PAIRING_FAILED 6
#define PARSE_FAILED 7
#define ALLOC_ERROR 8

/* Define data structures */
struct cmdparam
{
    int sa;
    int sb;
    int gap_open;
    int gap_extend;
    int min_score;
    int is_paired;
    int default_dir;
    char *filename1;
    char *filename2;
    char *parentdir;
    char *outdir;
};

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

char alpha[5] = "ACGTN";

/* Function prototypes */
struct cmdparam *parse_cmdline (int argc, char *argv[]);
void usage (char *program);

int
main (int argc, char *argv[])
{
    int i = 0;
    int j = 0;
    int k = 0;
    int xtra = KSW_XSTART;
    int writable = 0;
    int status = 0;
    unsigned int x = 0;
    int count = 0;
    char mat[25];
    struct cmdparam *cp;
    FQDB *reads;

    opterr = 0;

    /* Parse command line options */
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
            /* Test if output directory already exists and
             * whether it contains previous output files */
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

    /* Initialize the scoring matrix */
    for (i = k = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
                {
                    mat[k++] = (i == j) ? cp->sa : -cp->sb;
                }

            /* Ambiguous base */
            mat[k++] = 0;
        }

    for (j = 0; j < 5; j++)
        {
            mat[k++] = 0;
        }

    /* Read fastQ sequences from input file */
    if (cp->is_paired)
        {
            reads = read_fastq (2, cp->filename1, cp->filename2);
        }
    else
        {
            reads = read_fastq (1, cp->filename1);
        }
    if (reads == NULL)
        {
            fputs ("Error read fastQ input.\n", stderr);
            return DATABASE_CREATION_FAILED;
        }

	/* Pair the mates in the database */
	if (cp->is_paired)
		{
			char *filename_base = NULL;
			char *p = NULL;
			char *cc = NULL;
			char *reserve = NULL;
			char seps[] = "_.";

			/* Allocate memory for filename_base */
			if ((filename_base = malloc (strlen (cp->outdir) + strlen (cp->filename1) + 2)) == NULL)
				{
					fprintf (stderr, "[trim_3prime:%s:%d] Error allocating memory "
							 "for filename_base.\n", __func__, __LINE__);
					return ALLOC_ERROR;
				}

			/* Parse input file name to get the sample ID */
			cc = strdup (cp->filename1);
			p = strtok_r (cc, seps, &reserve);
			p = strtok_r (NULL, seps, &reserve);
			strcpy (filename_base, cp->outdir);
			strcat (filename_base, "final_");
			strcat (filename_base, p);
			free (cc);

			if ((pair_fastq_mates (reads)) < 0)
				{
					fprintf (stderr, "[trim_3prime:%s:%d] Error pairing "
					         "mates in fastQ database.\n", __func__,
					         __LINE__);
					return PAIRING_FAILED;
				}

			/* Trim 3' end of reverse reads that have adapter sequence */
			for (x = fqdb_begin(reads); x != fqdb_end(reads); x++)
				{
					if (fqdb_exists(reads, x) &&
					    fqdb_value(reads, x)->has_mate &&
					    fqdb_value(reads, x)->read == 1)
						{
							ALIGN_QUERY *q[2] = {0, 0};
							ALIGN_RESULT r;
							unsigned int m = fqdb_value(reads, x)->mate;
							char *target;
							char *query;
							target = strdup (fqdb_value(reads, x)->seq);
							query = revcom (fqdb_value(reads, m)->seq);
							size_t tlen = strlen (target);
							size_t qlen = strlen (query);

							/* Transform sequences */
							for (i = 0; i < qlen; i++)
								{
									query[i] = seq_nt4_table[(int)query[i]];
								}

							for (i = 0; i < tlen; i++)
								{
									target[i] = seq_nt4_table[(int)target[i]];
								}

							/* Do the alignment */
							r = local_align ((int)qlen,
							                 (unsigned char *)query,
										     (int)tlen,
										     (unsigned char *)target,
										     5, mat, cp->gap_open,
										     cp->gap_extend, xtra,
										     &q[0]);

							if (r.score >= cp->min_score)
								{
									/* Test trimming criterion */
									if ((r.target_begin == 0) &&
									    (r.query_begin > 0))
										{
											int new_end_pos = qlen - r.query_begin;
											char *seq = fqdb_value(reads, m)->seq;
											char *qual = fqdb_value(reads, m)->qual;
											seq[new_end_pos] = '\0';
											qual[new_end_pos] = '\0';
											seq = realloc (seq, (new_end_pos + 1) * sizeof (char));
											qual = realloc (qual, (new_end_pos + 1) * sizeof (char));
											count++;
										}
								}

							free (q[0]);
							free (q[1]);
							free (target);
							free (query);
						}
				}
			write_fastq_paired (reads, filename_base);
		}

	free (cp->filename1);
	if (cp->is_paired)
		{
			free (cp->filename2);
		}
	if (!cp->default_dir)
		{
			free (cp->parentdir);
		}
	free (cp->outdir);
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
			fprintf (stderr, "[trim_3prime:%s:%d] Error allocating memory "
					 "for command line structure.\n", __func__, __LINE__);
			return NULL;
		}

	/* Initialize default values in command line option structure */
	cp->sa = 1;
	cp->sb = 3;
	cp->gap_open = 5;
	cp->gap_extend = 2;
	cp->min_score = 100;
	cp->is_paired = 0;
	cp->filename1 = NULL;
	cp->filename2 = NULL;
	cp->default_dir = 0;
	cp->parentdir = NULL;
	cp->outdir = NULL;

	/* parse command line */
	while ((c = getopt (argc, argv, "a:b:q:o:hr:t:")) >= 0)
		{
			switch (c)
				{
					case 'a':
						cp->sa = atoi (optarg);
						break;
					case 'b':
						cp->sb = atoi (optarg);
						break;
					case 'o':
                        cp->parentdir = strdup (optarg);
                        size = strlen (cp->parentdir);
                        if (cp->parentdir[size - 1] == '/')
                            {
                                cp->outdir = malloc (size + 7);
                                strcpy (cp->outdir, cp->parentdir);
                                strcat (cp->outdir, "final/");
                            }
                        else
                            {
                                cp->outdir = malloc (size + 8);
                                strcpy (cp->outdir, cp->parentdir);
                                strcat (cp->outdir, "/final/");
                            }
						break;
					case 'q':
						cp->gap_open = atoi (optarg);
						break;
					case 'r':
						cp->gap_extend = atoi (optarg);
						break;
					case 't':
						cp->min_score = atoi (optarg);
						break;
					case 'h':
						return NULL;
					case '?':
						if (optopt == 'c')
							{
								fprintf (stderr, "Option -%c requires "
								         "an argument.\n", optopt);
							}
						else if (isprint (optopt))
							{
								fprintf (stderr, "Unknown option `-%c'.\n",
								         optopt);
							}
						else
							{
								fprintf (stderr, "Unknown option character "
								         "`\\x%x'.\n", optopt);
							}
						return NULL;
					default:
						return NULL;
				}
		}

	/* Get non-option arguments */
	if ((optind + 1) > argc)
		{
			return NULL;
		}
	else
		{
			cp->filename1 = strdup (argv[optind]);
			if ((optind + 2) <= argc)
				{
					cp->is_paired = 1;
					cp->filename2 = strdup (argv[optind + 1]);
				}
            if (cp->parentdir == NULL)
                {
                    char *fullpath = strdup (argv[optind]);
                    cp->parentdir = dirname (fullpath);
                    size = strlen (cp->parentdir);
					cp->outdir = malloc (size + 8);
					strcpy (cp->outdir, cp->parentdir);
					strcat (cp->outdir, "/final/");
					cp->default_dir = 1;
                }
			return cp;
		}
}

void
usage (char *program)
{
    fprintf (stderr, "Usage : %s [OPTIONS] [FASTQ.R1] [FASTQ.R2]\n", program);
    fputs ("Aligns mate pairs in fastQ files and trims 3' end of "
           "reverse sequence\n\n", stderr);
    fputs ("Available options\n", stderr);
    fputs ("  -a  INT    Log odds score for identical bases "
           "[default: 1]\n", stderr);
    fputs ("  -b  INT    Log odds score for different bases "
           "[default: 3]\n", stderr);
    fputs ("  -o  DIR    Parent directory to write output files "
           "[default: same as input fastQ]\n", stderr);
    fputs ("  -q  INT    Gap open cost [default: 5]\n", stderr);
    fputs ("  -r  INT    Gap extend cost [default: 2]\n", stderr);
    fputs ("  -t  INT    Minimum alignment score [default: 100]\n",
           stderr);
    fputs ("  -h         Display this help message\n\n", stderr);
}
