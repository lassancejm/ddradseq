/* file: print_buffer.c
 * description: Function to print current buffer to stdout
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

int
print_buffer(char *buff, const size_t nl)
{
	char *p = buff;
	size_t i = 0;
	size_t l = 0;

	for (i = 0; i < nl; i++)
	{
		l = strlen(p);
		fprintf(stdout, "%s\n", p);
		p += l + 1;
	}
	return 0;
}
