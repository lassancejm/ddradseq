/* file: reset_buffer.c
 * description: Function to reset the buffer for next read block
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include "ddradseq.h"

size_t
reset_buffer(char *buff, const char *r)
{
	size_t i = 0;
	size_t br = 0;

	br = strlen(r);
	for (i = 0; i < br; i++)
		buff[i] = *(r + i);
	return br;
}
