#include <stdio.h>
#include <stdlib.h>
#include <bzlib.h>

#define BUFLEN 16384

int main(int argc, char *argv[])
{
	FILE* f;
	BZFILE* b;
	int nBuf;
	char buf[BUFLEN];
	int bzerror;
	int nWritten;

	f = fopen("test.R1.fq.bz2", "r");
	if (!f)
	{
		fputs("Failed to open file.\n", stderr);
		return 1;
	}
	b = BZ2_bzReadOpen(&bzerror, f, 0, 0, NULL, 0);
	if (bzerror != BZ_OK)
	{
		BZ2_bzReadClose(&bzerror, b);
		fputs("bzip2 error.\n", stderr);
		return 1;
	}

	bzerror = BZ_OK;
	while (bzerror == BZ_OK)
	{
		nBuf = BZ2_bzRead(&bzerror, b, buf, BUFLEN);
		if (bzerror == BZ_OK || bzerror == BZ_STREAM_END)
		{
			int i;
			for (i=0; i < nBuf; i++) putchar(buf[i]);
		}
	}

	if (bzerror != BZ_STREAM_END)
	{
		BZ2_bzReadClose(&bzerror, b);
		fputs("bzip error after read.\n", stderr);
		return 1;
	}
	else
	{
		BZ2_bzReadClose(&bzerror, b);
	}
	fclose(f);
	return 0;
}
