/* file: read_csv.c
 * description: Function to read CSV file into hash database
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "khash.h"
#include "ddradseq.h"

#define MAX_LINE 400

khash_t(pool_hash) *
read_csv (CMD *cp)
{
	char *csvfile = cp->csvfile;  /* Pointer to CSV database file name */
	char *outpath = cp->outdir;	  /* Pointer to parent of output directories */
	char buf[MAX_LINE];			  /* File input buffer */
	char seps[] = ",";			  /* CSV entry separator character */
	char *tok = NULL;			  /* Holds parsed CSV tokens */
	char *r = NULL;				  /* Residual pointer for strtok_r */
	char *tmp = NULL;			  /* Temporary pointer */
	unsigned char trail = 0;	  /* Boolean indicator of trailing slash */
	int a = 0;					  /* Return value for database entry */
	size_t strl = 0;			  /* Generic string length holder */
	size_t pathl = 0;			  /* Length of path string */
	gzFile in;					  /* Input file stream */
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khash_t(barcode) *b = NULL;
	khash_t(pool) *p = NULL;
	khash_t(pool_hash) *h = NULL;
	BARCODE *bc = NULL;
	POOL *pl = NULL;

	/* Check for trailing slash on outpath */
	pathl = strlen(outpath);
	if (outpath[pathl - 1u] == '/') trail = 1;

	/* Open input database text file stream */
	if ((in = gzopen(csvfile, "rb")) == NULL)
	{
		fprintf(stderr, "Error opening database file \'%s\'.\n", csvfile);
		return NULL;
	}

	/* Initialize top-level hash */
	h = kh_init(pool_hash);

	/* Read CSV and populate the database */
	while (gzgets(in, buf, MAX_LINE) != Z_NULL)
	{
		/* Re-initialize measure of outfile path string length */
		pathl = strlen(outpath);

		/* Get the flowcell entry */
		if ((tok = strtok_r(buf, seps, &r)) == NULL)
		{
			fputs("ERROR: strtok_r failed.\n", stderr);
			exit(EXIT_FAILURE);
		}

		/* Alloc memory for flowcell string */
		strl = strlen(tok);
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for tmp.\n", stderr);
			exit(EXIT_FAILURE);
		}
		pathl += strl + 1u;
		strcpy(tmp, tok);

		/* Put flowcell string as key in top-level hash */
		i = kh_put(pool_hash, h, tmp, &a);

		/* If this flowcell is a new entry-- */
		/* initialize a second-level hash and add to value */
		if (a)
		{
			p = kh_init(pool);
			kh_value(h, i) = p;
		}
		else
			free(tmp);

		/* Get pool sequence */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			fputs("ERROR: strtok_r failed.\n", stderr);
			exit(EXIT_FAILURE);
		}

		/* Alloc memory for pool sequence string */
		strl = strlen(tok);
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for tmp.\n", stderr);
			exit(EXIT_FAILURE);
		}
		strcpy(tmp, tok);

		/* Put pool sequence string in second-level hash */
		p = kh_value(h, i);
		j = kh_put(pool, p, tmp, &a);

		/* If this pool sequence is a new entry-- */
		/* initialize a third-level hash and POOL data structure */
		/* Add the new hash to POOL and add POOL to second-level hash */
		if (a)
		{
			if ((pl = malloc(sizeof(POOL))) == NULL)
			{
				fputs("ERROR: cannot allocate memory for pl.\n", stderr);
				exit(EXIT_FAILURE);
			}
			b = kh_init(barcode);
			pl->b = b;
			kh_value(p, j) = pl;
		}
		else
			free(tmp);

		/* Get pool value */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			fputs("ERROR: strtok_r failed.\n", stderr);
			exit(EXIT_FAILURE);
		}

		/* Alloc memory for pool identifier */
		strl = strlen(tok);
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for tmp.\n", stderr);
			exit(EXIT_FAILURE);
		}
		pathl += strl + 1u;
		strcpy(tmp, tok);

		/* Put pool identifier string and directory path */
		/* in POOL data structure */
		pl = kh_value(p, j);
		pl->poolID = tmp;
		if ((tmp = malloc(pathl + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for tmp.\n", stderr);
			exit(EXIT_FAILURE);
		}
		if (trail)
			sprintf(tmp, "%s%s/%s", outpath, kh_key(h, i), pl->poolID);
		else
			sprintf(tmp, "%s/%s/%s", outpath, kh_key(h, i), pl->poolID);
		pl->poolpath = tmp;

		/* Get barcode sequence */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			fputs("ERROR: strtok_r failed.\n", stderr);
			exit(EXIT_FAILURE);
		}

		/* Alloc memory for barcode sequence string */
		strl = strlen(tok);
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for tmp.\n", stderr);
			exit(EXIT_FAILURE);
		}
		strcpy(tmp, tok);

		/* Put barcode sequence in third-level hash */
		pl = kh_value(p, j);
		b = pl->b;
		k = kh_put(barcode, b, tmp, &a);
		if (b->size == 1) pl->barcode_length = strl;
		else
		{
			if (pl->barcode_length != strl)
			{
				fputs("Error: unequal barcode lengths\n", stderr);
				abort();
			}
		}
		/* If this barcode is a new entry-- */
		/* initialize a fourth-level hash and add to value */
		if (a)
		{
			if ((bc = malloc(sizeof(BARCODE))) == NULL)
			{
				fputs("ERROR: cannot allocate memory for bc.\n", stderr);
				exit(EXIT_FAILURE);
			}
			if ((bc->buffer = malloc(BUFLEN)) == NULL)
			{
				fputs("ERROR: cannot allocate memory for bc->buffer.\n", stderr);
				exit(EXIT_FAILURE);
			}
			bc->buffer[0] = '\0';
			bc->curr_bytes = 0;
			kh_value(b, k) = bc;
		}
		else
			free(tmp);

		/* Get barcode value */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			fputs("ERROR: strtok_r failed.\n", stderr);
			exit(EXIT_FAILURE);
		}

		/* Alloc memory for barcode value */
		strl = strcspn(tok, " \n");
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for tmp.\n", stderr);
			exit(EXIT_FAILURE);
		}
		pathl += strl + 1u;
		strncpy(tmp, tok, strl);
		tmp[strl] = '\0';

		/* Add barcode value to BARCODE data structure */
		bc = kh_value(b, k);
		bc->smplID = tmp;
		pathl += 21u;
		if ((tmp = malloc(pathl + 1u)) == NULL)
		{
			fputs("ERROR: cannot allocate memory for tmp.\n", stderr);
			exit(EXIT_FAILURE);
		}
		sprintf(tmp, "%s/parse/smpl_%s.R1.fq.gz", pl->poolpath, bc->smplID);
		bc->outfile = tmp;
	}

	/* Close input CSV file stream */
	gzclose(in);

	return h;
}
