/* file: traverse_dirtree.c
 * description: Function to get a sorted list of all fastQ files in a directory tree
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
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

#ifndef USE_FDS
#define USE_FDS 15
#endif

/* Globally scoped variables */
unsigned int n;
char **f;

/* Function prototypes */
int count_parsefiles(const char *filepath, const struct stat *info,
					 const int typeflag, struct FTW *pathinfo);
int get_parsefiles(const char *filepath, const struct stat *info,
				   const int typeflag, struct FTW *pathinfo);
int count_pairfiles(const char *filepath, const struct stat *info,
					const int typeflag, struct FTW *pathinfo);
int get_pairfiles(const char *filepath, const struct stat *info,
				  const int typeflag, struct FTW *pathinfo);
static int compare(const void * a, const void * b);

char**
traverse_dirtree(const char *dirpath, char *pattern, unsigned int *x)
{
	int r = 0;
	unsigned int i = 0;

	/* Check validity of directory path */
	if (dirpath == NULL || *dirpath == '\0')
		return NULL;

	/* Count number of files in directory tree */
	n = 0;
	if (strcmp(pattern, "parse") == 0)
		r = nftw(dirpath, count_parsefiles, USE_FDS, FTW_PHYS);
	else if (strcmp(pattern, "pairs") == 0)
		r = nftw(dirpath, count_pairfiles, USE_FDS, FTW_PHYS);

	/* Check for errors */
	if (r >= 0)
		errno = r;
	if (errno)
	{
		fprintf(stderr, "%s.\n", strerror(errno));
		return NULL;
	}

	/* Get list of filenames */
	f = malloc(n * sizeof(char*));
	n = 0;
	if (strcmp(pattern, "parse") == 0)
		r = nftw(dirpath, get_parsefiles, USE_FDS, FTW_PHYS);
	else if (strcmp(pattern, "pairs") == 0)
		r = nftw(dirpath, get_pairfiles, USE_FDS, FTW_PHYS);

	/* Check for errors */
	if (r >= 0)
		errno = r;
	if (errno)
	{
		fprintf(stderr, "%s.\n", strerror(errno));
		return NULL;
	}

	/* Assign number of files */
	*x = n;

	/* Sort file list */
	qsort(f, n, sizeof(const char *), compare);

	return f;
}

int
count_parsefiles(const char *filepath, const struct stat *info, const int typeflag,
				 struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "parse");
		if (p != NULL && q != NULL) n++;
	}
	return 0;
}

int
get_parsefiles(const char *filepath, const struct stat *info, const int typeflag,
			   struct FTW *pathinfo)
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

int
count_pairfiles(const char *filepath, const struct stat *info, const int typeflag,
				struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "pairs");
		if (p != NULL && q != NULL) n++;
	}
	return 0;
}

int
get_pairfiles(const char *filepath, const struct stat *info, const int typeflag,
			  struct FTW *pathinfo)
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

static int
compare(const void * a, const void * b)
{
	return strcmp(*(const char **) a, *(const char **) b);
}
