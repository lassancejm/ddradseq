/* file: error.c
 * description: Error reporting function for ddradseq
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdarg.h>
#include "ddradseq.h"

void logerror(const char *format, ...)
{
	va_list ap;
	va_list copy;

	/* Update time string */
	get_timestr(&timestr[0]);

	fprintf(stderr, "[ddradseq: %s] ERROR -- ", timestr);
	fprintf(lf, "[ddradseq: %s] ERROR -- ", timestr);

	va_start(ap, format);
	va_copy(copy, ap);
	vfprintf(stderr, format, ap);
	vfprintf(lf, format, copy);
	va_end(ap);
}

void error(const char *format, ...)
{
	va_list ap;

	/* Update time string */
	get_timestr(&timestr[0]);

	fprintf(stderr, "[ddradseq: %s] ERROR -- ", timestr);

	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
}
