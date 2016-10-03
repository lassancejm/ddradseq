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
	khint_t i = 0;                /* Generic hash iterator */
	khint_t j = 0;                /* Generic hash iterator */
	khint_t k = 0;                /* Generic hash iterator */
	khash_t(barcode) *b = NULL;   /* Pointer to barcode hash table */
	khash_t(pool) *p = NULL;      /* Pointer to pool hash table */
	khash_t(pool_hash) *h = NULL; /* Pointer to flow hash table */
	BARCODE *bc = NULL;           /* Pointer barcode data structure */
	POOL *pl = NULL;              /* Pointer to pool data structure */

	/* Get time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log */
	fprintf(lf, "[ddradseq: %s] INFO -- Parsing CSV database file \'%s\'.\n",
	        timestr, csvfile);

	/* Check for trailing slash on outpath */
	pathl = strlen(outpath);
	if (outpath[pathl - 1u] == '/') trail = 1;

	/* Open input database text file stream */
	if ((in = gzopen(csvfile, "rb")) == NULL)
	{
		fprintf(stderr, "ERROR: Could not open database file \'%s\'.\n", csvfile);
		fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Could not open CSV database file "
		        " %s into memory.\n", timestr, __func__, __LINE__, csvfile);
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
			fputs("ERROR: Parsing CSV failed.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing CSV file failed.\n",
			        timestr, __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for flowcell string */
		strl = strlen(tok);
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			return NULL;
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
			fputs("ERROR: Parsing CSV failed.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing CSV file failed.\n",
			        timestr, __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for pool sequence string */
		strl = strlen(tok);
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			return NULL;
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
				fputs("ERROR: Memory allocation failure.\n", stderr);
				fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
				        timestr, __func__, __LINE__);
				return NULL;
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
			fputs("ERROR: Parsing CSV failed.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing CSV file failed.\n",
			        timestr, __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for pool identifier */
		strl = strlen(tok);
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			return NULL;
		}
		pathl += strl + 1u;
		strcpy(tmp, tok);

		/* Put pool identifier string and directory path */
		/* in POOL data structure */
		pl = kh_value(p, j);
		pl->poolID = tmp;
		if ((tmp = malloc(pathl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			return NULL;
		}
		if (trail)
			sprintf(tmp, "%s%s/%s", outpath, kh_key(h, i), pl->poolID);
		else
			sprintf(tmp, "%s/%s/%s", outpath, kh_key(h, i), pl->poolID);
		pl->poolpath = tmp;

		/* Get barcode sequence */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			fputs("ERROR: Parsing CSV failed.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing CSV file failed.\n",
			        timestr, __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for barcode sequence string */
		strl = strlen(tok);
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			return NULL;
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
				fputs("ERROR: Unequal barcode lengths.\n", stderr);
				fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Unequal barcode lengths "
				        "in CSV file %s.\n", timestr, __func__, __LINE__, csvfile);
				return NULL;
			}
		}
		/* If this barcode is a new entry-- */
		/* initialize a fourth-level hash and add to value */
		if (a)
		{
			if ((bc = malloc(sizeof(BARCODE))) == NULL)
			{
				fputs("ERROR: Memory allocation failure.\n", stderr);
				fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
				        timestr, __func__, __LINE__);
				return NULL;
			}
			if ((bc->buffer = malloc(BUFLEN)) == NULL)
			{
				fputs("ERROR: Memory allocation failure.\n", stderr);
				fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
				        timestr, __func__, __LINE__);
				return NULL;
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
			fputs("ERROR: Parsing CSV failed.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Parsing CSV file failed.\n",
			        timestr, __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for barcode value */
		strl = strcspn(tok, " \n");
		if ((tmp = malloc(strl + 1u)) == NULL)
		{
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			return NULL;
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
			fputs("ERROR: Memory allocation failure.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Memory allocation failure.\n",
			        timestr, __func__, __LINE__);
			return NULL;
		}
		sprintf(tmp, "%s/parse/smpl_%s.R1.fq.gz", pl->poolpath, bc->smplID);
		bc->outfile = tmp;
	}

	/* Close input CSV file stream */
	gzclose(in);

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Print informational message to log */
	fprintf(lf, "[ddradseq: %s] INFO -- Successfully parsed CSV database file \'%s\'.\n",
	        timestr, csvfile);

	return h;
}
