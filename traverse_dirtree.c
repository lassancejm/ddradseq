/* file: traverse_dirtree.c
 * description: Produces a sorted list of all fastQ files in the input directory tree
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <ftw.h>
#include "ddradseq.h"

#ifndef USE_FDS
#define USE_FDS 15
#endif

/* Globally scoped variables */
unsigned int n;
char **f;
extern int errno;

/* Function prototypes */
static int count_fastqfiles(const char *filepath, const struct stat *info,
                            const int typeflag, struct FTW *pathinfo);
static int get_fastqfiles(const char *filepath, const struct stat *info,
                          const int typeflag, struct FTW *pathinfo);
static int count_parsefiles(const char *filepath, const struct stat *info,
	                        const int typeflag, struct FTW *pathinfo);
static int get_parsefiles(const char *filepath, const struct stat *info,
                          const int typeflag, struct FTW *pathinfo);
static int count_pairfiles(const char *filepath, const struct stat *info,
                           const int typeflag, struct FTW *pathinfo);
static int get_pairfiles(const char *filepath, const struct stat *info,
                         const int typeflag, struct FTW *pathinfo);
static int compare(const void * a, const void * b);

unsigned int traverse_dirtree(const CMD *cp, char **flist)
{
	char *errstr = NULL;
	const char *dirpath = (string_equal(cp->mode, "pair") || string_equal(cp->mode, "all")) ? cp->parent_indir : cp->outdir;
	int r = 0;
	unsigned int x = 0;
	FILE *lf = cp->lf;

	/* Check validity of directory path */
	if (dirpath == NULL || *dirpath == '\0')
		return 0;

	/* Count number of files in directory tree */
	n = 0;
	if (string_equal(cp->mode, "pair"))
		r = nftw(dirpath, count_parsefiles, USE_FDS, FTW_PHYS);
	else if (string_equal(cp->mode, "trimend"))
		r = nftw(dirpath, count_pairfiles, USE_FDS, FTW_PHYS);
	else
		r = nftw(dirpath, count_fastqfiles, USE_FDS, FTW_PHYS);

	/* Check for errors */
	if (r < 0)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Directory traversal on %s failed: %s.\n", __func__, __LINE__,
		         dirpath, errstr);
		return 0;
	}
	else if (r > 0)
	{
		logerror(lf, "%s:%d Directory traversal callback failed.\n", __func__, __LINE__);
		return 0;
	}

	/* Get list of filenames */
	f = malloc(n * sizeof(char*));
	if (UNLIKELY(!f))
	{
		logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 0;
	}
	n = 0;
	if (string_equal(cp->mode, "pair"))
		r = nftw(dirpath, get_parsefiles, USE_FDS, FTW_PHYS);
	else if (string_equal(cp->mode, "trimend"))
		r = nftw(dirpath, get_pairfiles, USE_FDS, FTW_PHYS);
	else
		r = nftw(dirpath, get_fastqfiles, USE_FDS, FTW_PHYS);

	/* Check for errors */
	if (r < 0)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Directory traversal failed: %s.\n", __func__, __LINE__,
		         errstr);
		return 0;
	}
	else if (r > 0)
	{
		logerror(lf, "%s:%d Directory traversal callback failed.\n", __func__, __LINE__);
		return 0;
	}

	/* Sort file list */
	qsort(f, n, sizeof(const char *), compare);

	/* Assign number of files */
	x = n;

	/* Assign address of file list */
	flist = f;

	return x;
}

static int count_fastqfiles(const char *filepath, const struct stat *info,
                            const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq");
		q = strstr(filepath, ".fastq");
		if (p != NULL || q != NULL)
			n++;
	}
	return 0;
}

static int get_fastqfiles(const char *filepath, const struct stat *info,
                          const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;
	size_t l = 0;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq");
		q = strstr(filepath, ".fastq");
		if (p != NULL || q != NULL)
		{
			l = strlen(filepath);
			f[n] = malloc(l + 1u);
			if (UNLIKELY(!f[n]))
			{
				error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
				return 1;
			}
			strcpy(f[n], filepath);
			n++;
		}
	}
	return 0;
}


static int count_parsefiles(const char *filepath, const struct stat *info,
                            const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "parse");
		if (p != NULL && q != NULL)
			n++;
	}
	return 0;
}

static int get_parsefiles(const char *filepath, const struct stat *info,
                          const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;
	size_t l = 0;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "parse");
		if (p != NULL && q != NULL)
		{
			l = strlen(filepath);
			f[n] = malloc(l + 1u);
			if (UNLIKELY(!f[n]))
			{
				error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
				return 1;
			}
			strcpy(f[n], filepath);
			n++;
		}
	}
	return 0;
}

static int count_pairfiles(const char *filepath, const struct stat *info,
                           const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "pairs");
		if (p != NULL && q != NULL)
			n++;
	}
	return 0;
}

static int get_pairfiles(const char *filepath, const struct stat *info,
                         const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;
	size_t l = 0;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "pairs");
		if (p != NULL && q != NULL)
		{
			l = strlen(filepath);
			f[n] = malloc(l + 1u);
			if (UNLIKELY(!f[n]))
			{
				error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
				return 1;
			}
			strcpy(f[n], filepath);
			n++;
		}
	}
	return 0;
}

static int compare(const void * a, const void * b)
{
	return strcmp(*(const char **) a, *(const char **) b);
}
