/* file: traverse_dirtree.c
 * description: Produces a sorted list of all fastQ files in the input directory tree
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
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

char **traverse_dirtree(const char *dirpath, const char *pattern, unsigned int *x)
{
	char *p = (char*)pattern;
	int r = 0;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Check validity of directory path */
	if (dirpath == NULL || *dirpath == '\0')
		return NULL;

	/* Count number of files in directory tree */
	n = 0;
	if (!p)
		r = nftw(dirpath, count_fastqfiles, USE_FDS, FTW_PHYS);
	else if (string_equal(p, "parse"))
		r = nftw(dirpath, count_parsefiles, USE_FDS, FTW_PHYS);
	else if (string_equal(p, "pairs"))
		r = nftw(dirpath, count_pairfiles, USE_FDS, FTW_PHYS);

	/* Check for errors */
	if (r >= 0)
		errno = r;
	if (errno)
	{
		logerror("%s:%d Directory traversal failed.\n", __func__, __LINE__);
		return NULL;
	}

	/* Get list of filenames */
	f = malloc(n * sizeof(char*));
	if (UNLIKELY(!f))
	{
		logerror("%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return NULL;
	}
	n = 0;
	if (!p)
		r = nftw(dirpath, get_fastqfiles, USE_FDS, FTW_PHYS);
	else if (string_equal(p, "parse"))
		r = nftw(dirpath, get_parsefiles, USE_FDS, FTW_PHYS);
	else if (string_equal(p, "pairs"))
		r = nftw(dirpath, get_pairfiles, USE_FDS, FTW_PHYS);

	/* Check for errors */
	if (r >= 0)
		errno = r;
	if (errno)
	{
		logerror("%s:%d Directory traversal failed.\n", __func__, __LINE__);
		return NULL;
	}

	/* Assign number of files */
	*x = n;

	/* Sort file list */
	qsort(f, n, sizeof(const char *), compare);

	return f;
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
