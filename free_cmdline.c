/* file: free_cmdline.c
 * description: Deallocate memory for command line parameter data structure
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include <stdbool.h>
#include "ddradseq.h"

int
free_cmdline(CMD *cp)
{
	if (cp->default_dir == false)
		free(cp->parentdir);
	free(cp->outdir);
	free(cp->filename);
	free(cp->csvfile);
	free(cp);
	return 0;
}