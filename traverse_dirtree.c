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

unsigned int n;
char **f;

int count_files(const char *filepath, const struct stat *info,
                const int typeflag, struct FTW *pathinfo);
int get_files(const char *filepath, const struct stat *info,
              const int typeflag, struct FTW *pathinfo);

char**
traverse_dirtree(const char *dirpath, unsigned int *x)
{
    int r = 0;
    unsigned int i = 0;

    /* Check validity of directory path */
    if (dirpath == NULL || *dirpath == '\0')
        return NULL;

    /* Count number of files in directory tree */
    n = 0;
    r = nftw(dirpath, count_files, USE_FDS, FTW_PHYS);
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
    r = nftw(dirpath, get_files, USE_FDS, FTW_PHYS);
    if (r >= 0)
        errno = r;
    if (errno)
    {
        fprintf(stderr, "%s.\n", strerror(errno));
        return NULL;
    }
    *x = n;

    return f;
}

int count_files(const char *filepath, const struct stat *info,
                const int typeflag, struct FTW *pathinfo)
{
    if (typeflag == FTW_F) n++;

    return 0;
}

int get_files(const char *filepath, const struct stat *info,
                const int typeflag, struct FTW *pathinfo)
{
    size_t l = 0;

    if (typeflag == FTW_F)
    {
        l = strlen(filepath);
        f[n] = malloc(l + 1u);
        strcpy(f[n], filepath);
        n++;
    }
    return 0;
}
