/* file: clean_buffer.c
 * description: Function to limit buffer to hold only entire fastQ entries
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

char *
clean_buffer(char *buff, size_t *nl)
{
	char *p = NULL;
	char *s = buff;
	size_t i = 0;

	p = strpbrk(s, "\n");
	if (p == NULL)
		return NULL;
	while (p)
	{
		if (i == (*nl - (*nl % 4)))
			break;
		*p = '\0';
		s = p + 1;
		p = strpbrk(s, "\n");
		i++;
	}
	*nl = i;
	return s;
}
