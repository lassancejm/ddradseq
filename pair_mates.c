/* file pair_mates.c
 * brief Pair mates in two fastQ files
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "ddradseq.h"
#include "khash.h"

#define MAX_LINE_LENGTH 400
#define BSIZE 4000
#define KEYLEN 31

int
pair_mates(char *filename, khash_t(fastq) *h, char *ffor, char *frev)
{
    char buf[BSIZE][MAX_LINE_LENGTH];
	char *idline = NULL;
    char seps[] = ": ";
    char *tok = NULL;
    char *mkey = NULL;
    char *r = NULL;
    int tile = 0;
    int xpos = 0;
    int ypos = 0;
    size_t z = 0;
    size_t l = 0;
    size_t lc = 0;
    size_t pos = 0;
    size_t strl = 0;
    khint_t k = 0;
    gzFile in;
    gzFile fout;
    gzFile rout;
    FASTQ *e = NULL;

    /* Open the fastQ input stream */
    if ((in = gzopen(filename, "rb")) == Z_NULL)
    {
        fprintf(stderr, "Error: cannot open input fastQ file %s.\n", filename);
        return 1;
    }

    /* Open the output fastQ file streams */
     if ((fout = gzopen(ffor, "wb")) == Z_NULL)
    {
        fprintf(stderr, "Error: cannot open input fastQ file %s.\n", ffor);
        return 1;
    }

    if ((rout = gzopen(frev, "wb")) == Z_NULL)
    {
        fprintf(stderr, "Error: cannot open input fastQ file %s.\n", frev);
        return 1;
    }

    /* Enter data from the fastQ input file into the database */
    while (1)
    {
        /* Fill up the buffer */
        for (lc = 0; lc < BSIZE; lc++)
        {
            /* Clear the buffer */
            memset(buf[lc], 0, MAX_LINE_LENGTH);

            /* Get line from the fastQ input stream */
            if (gzgets(in, buf[lc], MAX_LINE_LENGTH) == Z_NULL) break;
        }

        /* Iterate through lines in the buffer */
        for (l = 0; l < lc; l++)
        {
            /* We are at the end of one fastQ entry */
            if (l % 4 == 3)
            {
                /* Parse entry identifier */
                pos = strcspn(buf[l - 3], "\n");
                buf[l - 3][pos] = '\0';
                strl = strlen(&buf[l - 3][1]);

                /* Construct fastQ hash key */
                idline = malloc(strl + 1u);
                assert(idline != NULL);
                strcpy(idline, &buf[l - 3][1]);

				/* Get instrument name */
				tok = strtok_r(idline, seps, &r);
				assert(tok != NULL);

				/* Get run number, flow cell ID, and lane number */
                for (z = 0; z < 3; z++)
                {
    				tok = strtok_r(NULL, seps, &r);
    				assert(tok != NULL);
                }

				/* Get tile number, x, and y coordinate */
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				tile = atoi(tok);
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				xpos = atoi(tok);
				tok = strtok_r(NULL, seps, &r);
				assert(tok != NULL);
				ypos = atoi(tok);
				mkey = malloc(KEYLEN);
				assert(mkey != NULL);
				sprintf(mkey, "%010d%010d%010d", tile, xpos, ypos);
				k = kh_get(fastq, h, mkey);
                if (k != kh_end(h))
                    e = kh_value(h, k);
                free(mkey);
				free(idline);

                if (e != NULL)
                {
                    /* Parse DNA sequence */
                    pos = strcspn(buf[l - 2], "\n");
                    buf[l - 2][pos] = '\0';
    
                    /* Parse quality sequence */
                    pos = strcspn(buf[l], "\n");
                    buf[l][pos] = '\0';

                    /* Need to construct output file streams */
                    gzprintf(fout, "@%s\n%s\n+\n%s\n", e->id, e->seq, e->qual);
                    gzprintf(rout, "@%s\n%s\n+\n%s\n", &buf[l - 3][1], &buf[l - 2][0], &buf[l][0]);
                }
            }
        }

        /* If we are at the end of the file */
        if (lc < BSIZE) break;
    }

    /* Close input stream */
    gzclose(in);
    gzclose(fout);
    gzclose(rout);

    return 0;
}