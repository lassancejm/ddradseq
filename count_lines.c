/* file: count_lines.c
 * description: Function to count newline characters in buffer
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

size_t
count_lines(const char *buff)
{
	char *p = NULL;
	size_t nl = 0;

	p = strchr(buff, '\n');
	while (p)
	{
		p = strchr(p + 1, '\n');
		nl++;
	}
	return nl;
}
