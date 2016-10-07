/* file: free_cmdline.c
 * description: Deallocates memory for command line parameter data structure
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include "ddradseq.h"

int
free_cmdline(CMD *cp)
{
	free(cp->parentdir);
	free(cp->outdir);
	free(cp->mode);
	free(cp->csvfile);
	free(cp);
	return 0;
}
