/* file: trim_barcode.c
 * description: Worker program for the ddRadSeq pipeline to trim custom
 * barcodes from 5' end of fastQ sequences
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <libgen.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pthread.h>
#include <zlib.h>
#include <fasta.h>

/* Define some main return errors */
#define INCORRECT_PARAMETERS 2
#define DIR_NOT_WRITABLE 3
#define CREATE_DIR_FAILED 4
#define DATABASE_CREATION_FAILED 5
#define ADAPTER_DB_FAILED 6
#define PAIRING_FAILED 7
#define ALLOC_ERROR 8

#define FORWARD 1
#define REVERSE 2

#define FALSE 0
#define TRUE 1

/* Define data structures */
struct cmdparam
{
    int dist;
    int nthreads;
    int is_paired;
    int default_dir;
    char *parentdir;
    char *outdir;
    char *filename1;
    char *filename2;
    char *csvfile;
};

struct output_file
{
    char *forname;
    char *revname;
    char *seq;
    gzFile ffp;
    gzFile rfp;
};

struct thread_data
{
    FQDB *readdata;
    struct output_file *outinfo;;
    int start;
    int end;
};

/* Function prototypes */
struct cmdparam *parse_cmdline (int argc, char *argv[]);
void usage (char *program);
void *pair_reads(void* data);

int
main (int argc, char *argv[])
{
    int writable = 0;
    int status = 0;
    unsigned int k = 0;
    struct cmdparam *cp;
    FQDB *reads;
    STDDB *adapt;

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
            fprintf (stderr, "[trim_barcode:%s:%d] Error: cannot write to "
                     "directory %s\n", __func__, __LINE__, cp->parentdir);
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
                            fprintf (stderr, "[trim_barcode[%s:%d] Error: "
                                     "failed to create output directory %s\n",
                                     __func__, __LINE__, cp->outdir);
                            return CREATE_DIR_FAILED;
                        }
                }
        }

    /* Read all fastQ sequences from input file into single database */
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
            fprintf (stderr, "[trim_barcode:%s:%d] Error read fastQ input\n",
                     __func__, __LINE__);
            return DATABASE_CREATION_FAILED;
        }

    /* Read the adapter CSV input file and create database */
    adapt = create_adapter_database (cp->csvfile);
    if (adapt == NULL)
        {
            fprintf (stderr, "[trim_barcode:%s:%d] Error cannot read CSV file "
                     "%s\n", __func__, __LINE__, cp->csvfile);
            return ADAPTER_DB_FAILED;
        }

    /* Pair all of the mates in the reads database */
    if (cp->is_paired)
        {
            if ((pair_fastq_mates (reads)) < 0)
                {
                    fprintf (stderr, "[trim_barcode:%s:%d] Error pairing mates "
                             "in fastQ database\n", __func__, __LINE__);
                    return PAIRING_FAILED;
                }

            /* Trim barcode sequence from 5' end of forward reads */
            trim_adapter_seq (reads, adapt, FORWARD, cp->dist, cp->nthreads);

            unsigned int adapt_size = hash_size(adapt);
            struct output_file *of = malloc (adapt_size * sizeof (struct output_file));
            int q = 0;

            /* Iterate through barcode sequences */
            for (k = hash_begin(adapt); k != hash_end(adapt); k++)
                {
                    if (hash_exists(adapt, k) && (hash_value(adapt, k) != NULL))
                        {
                            /* Assign outfile names */
                            of[q].seq = hash_value(adapt, k);
                            if ((of[q].forname = malloc (strlen (cp->outdir) + strlen (of[q].seq) + 16)) == NULL)
                                {
                                    fprintf (stderr, "[trim_barcode:%s:%d] Error: cannot allocate "
                                             "memory for outfile_for name.\n", __func__, __LINE__);
                                    return ALLOC_ERROR;
                                }
                            sprintf (of[q].forname, "%strim_%s.R1.fq.gz", cp->outdir, of[q].seq);
                            if ((of[q].revname = malloc (strlen (cp->outdir) + strlen (of[q].seq) + 16)) == NULL)
                                {
                                    fprintf (stderr, "[trim_barcode:%s:%d] Error: cannot allocate "
                                             "memory for outfile_rev name.\n", __func__, __LINE__);
                                    return ALLOC_ERROR;
                                }
                            sprintf (of[q].revname, "%strim_%s.R2.fq.gz", cp->outdir, of[q].seq);
                            of[q].ffp = gzopen (of[q].forname, "ab");
                            of[q].rfp = gzopen (of[q].revname, "ab");
                            q++;
                        }
                }

            /* Do multithreaded pairing of reads */
            int tasks_per_thread = (q + cp->nthreads - 1) / cp->nthreads;
            int x = 0;
            struct thread_data *data;
            pthread_t *threads;
            data = calloc (cp->nthreads, sizeof (struct thread_data));
            threads = malloc (cp->nthreads * sizeof (pthread_t));
            for (x = 0; x < cp->nthreads; x++)
                {
                    data[x].readdata = reads;
                    data[x].outinfo = &of[0];
                    data[x].start = x * tasks_per_thread;
                    data[x].end = (x + 1) * tasks_per_thread;
                }
            data[cp->nthreads - 1].end = q;
            for (x = 0; x < cp->nthreads; x++)
                {
                    pthread_create (threads+x, NULL, pair_reads, (void*)(&data[x]));
                }
            for (x = 0; x < cp->nthreads; x++)
                {
                    pthread_join (threads[x], NULL);
                }
            free (of);
        }   /* Done pairing mates in database */

    /* Deallocate memory */
    free (cp->filename1);
    free (cp->filename2);
    if (!cp->default_dir)
        {
            free (cp->parentdir);
        }
    free (cp->outdir);
    free (cp->csvfile);
    free (cp);

    return EXIT_SUCCESS;
}

void *pair_reads(void *data)
{
    int k;
    unsigned int j, m;
    struct thread_data *d = (struct thread_data*)data;
    int start = d->start;
    int end = d->end;
    FQDB *h = d->readdata;
    struct output_file *f = d->outinfo;

    for (k = start; k < end; k++)
        {
            /* Iterate through fastQ database */
            for (j = fqdb_begin(h); j != fqdb_end(h); j++)
                {
                    if (fqdb_exists(h, j) &&
                        (fqdb_value(h, j)->indiv_id != NULL) &&
                        (fqdb_value(h, j)->has_mate))
                        {
                            m = fqdb_value(h, j)->mate;
                            if ((strcmp (fqdb_value(h, j)->indiv_id, f[k].seq) == 0) &&
                                fqdb_exists(h, m))
                                {
                                    append_fastq_file (fqdb_value(h, j), fqdb_key(h, j), &f[k].ffp);
                                    append_fastq_file (fqdb_value(h, m), fqdb_key(h, m), &f[k].rfp);
                                }
                        }
                }
           gzclose (f[k].ffp);
           gzclose (f[k].rfp);
           free (f[k].revname);
           free (f[k].forname);
       }
}

struct cmdparam *
parse_cmdline (int argc, char *argv[])
{
    int c = 0;
    size_t size = 0;
    struct cmdparam *cp = NULL;

    if ((cp = malloc (sizeof (struct cmdparam))) == NULL)
        {
            fprintf (stderr, "[trim_barcode:%s:%d] Error: cannot allocate "
                     "memory for command line data structure\n",
                     __func__, __LINE__);
            return NULL;
        }

     /* Initialize to default values */
     cp->default_dir = FALSE;
     cp->parentdir = NULL;
     cp->outdir = NULL;
     cp->filename1 = NULL;
     cp->filename2 = NULL;
     cp->csvfile = NULL;
     cp->nthreads = 1;
     cp->is_paired = FALSE;
     cp->dist = 1;

    /* Parse command line arguments */
    while ((c = getopt (argc, argv, "o:d:t:h")) != -1)
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
                                strcat (cp->outdir, "trim/");
                            }
                        else
                            {
                                cp->outdir = malloc (size + 7);
                                strcpy (cp->outdir, cp->parentdir);
                                strcat (cp->outdir, "/trim/");
                            }
                        break;
                    case 't':
                        cp->nthreads = atoi (optarg);
                        break;
                    case 'd':
                        cp->dist = atoi (optarg);
                        break;
                    case 'h':
                        free (cp);
                        return NULL;
                    case '?':
                        if (optopt == 'c')
                            {
                                fprintf (stderr, "Option -%c requires "
                                         "an argument\n", optopt);
                            }
                        else if (isprint (optopt))
                            {
                                fprintf (stderr, "Unknown option "
                                         "`-%c'.\n", optopt);
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

    if ((optind + 3) > argc)
        {
            free (cp);
            return NULL;
        }
    else
        {
            cp->is_paired = TRUE;
            cp->filename1 = strdup (argv[optind]);
            cp->filename2 = strdup (argv[optind + 1]);
            cp->csvfile = strdup (argv[optind + 2]);
            if (cp->parentdir == NULL)
                {
                    char *fullpath = realpath (cp->filename1, NULL);
                    cp->parentdir = dirname (fullpath);
                    size = strlen (cp->parentdir);
                    cp->outdir = malloc (size + 7);
                    strcpy (cp->outdir, cp->parentdir);
                    strcat (cp->outdir, "/trim/");
                    cp->default_dir = TRUE;
                }
            return cp;
        }
}

void
usage (char *program)
{
    fprintf (stderr, "Usage : %s [OPTIONS] [FASTQ.R1] [FASTQ.R2] [CSV] \n", program);
    fputs ("Trim custom barcode sequences from fastQ file\n\n", stderr);
    fputs ("Available options\n", stderr);
    fputs ("  -d  INT    Edit distance [default: 1]\n", stderr);
    fputs ("  -o  DIR    Parent directory to write output files [default: same as input fastQ]\n", stderr);
    fputs ("  -t  INT    Number of threads to run [default: 1]\n", stderr);
    fputs ("  -h         Display this help message\n\n", stderr);
}
