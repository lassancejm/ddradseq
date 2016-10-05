/* file: revcom.c
 * description: Functions to reverse complement a DNA string with full IUPAC alphabet
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "ddradseq.h"

/* Beginning of ASCII DNA sequence alphabet */
#define DNA_BEGIN 65

/* Function prototypes */
int reverse_string(char *s);
int complement_string(char *s);

char *
revcom(char *input_string)
{
	char *str = NULL;
	char *pch = NULL;
	size_t strl = 0;
	size_t i = 0;

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Create local copy of DNA string */
	str = strdup(input_string);
	pch = strchr(str,'\n');
	if (pch) *pch = '\0';
	strl = strlen(str);

	/* Convert to upper-case characters */
	for (i = 0; i < strl; i++)
	{
		/* Check that characters are either alphabetical or gaps */
		if (isalpha(str[i]))
			str[i] = toupper(str[i]);
		else if (str[i] == '-')
			continue;
		else
		{
			fputs("ERROR: Reverse complementing sequence failed.\n", stderr);
			fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Bad character '\%c\' at position %zu.\n",
					timestr, __func__, __LINE__, str[i], i + 1u);
			return NULL;
		}
	}

	/* If string is not null and longer than one character, */
	/* then do the reverse complement */
	if (str && strlen(str) > 1u)
	{
		/* Reverse the string */
		reverse_string (str);

		/* Complement string */
		complement_string (str);
	}
	return str;
}

int
reverse_string(char *s)
{
	char *p1 = s;
	char *p2 = s + strlen(s) - 1;

	/* Reverse string */
	while (p2 > p1)
	{
		*p1 ^= *p2;
		*p2 ^= *p1;
		*p1 ^= *p2;
		++p1;
		--p2;
	}
	return 0;
}

int
complement_string(char *s)
{
	char j = 0;
	char *iupac = "ACGTURYSWKMN-";
	char *iupac_extend = "BDHV";
	size_t i = 0;
	unsigned int lookup_table[25] = { 6u, 0u, 19u,	0u, 0u, 0u, 0u,	 0u,  0u, 0u,
									 12u, 0u, 10u, 13u, 0u, 0u, 0u, 17u, 22u, 2u,
									  2u, 0u, 18u,	0u, 24u};

	/* Update time string */
	get_timestr(&timestr[0]);

	/* Iterate through string and complement each base */
	for (i = 0; i < strlen(s); i++)
	{
		j = s[i];

		/* Check that base is IUPAC representation */
		if (strchr(iupac, j) == NULL)
		{
			if (strchr(iupac_extend, j) != NULL)
			{
				fputs("ERROR: Base complementing failed.\n", stderr);
				fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d IUPAC codes with three "
				        "bases at a site are not supported.\n", timestr, __func__, __LINE__);
			}
			else
			{
				fputs("ERROR: Base complementing failed.\n", stderr);
				fprintf(lf, "[ddradseq: %s] ERROR -- %s:%d Bad character '\%c\' at position %zu.\n",
					    timestr, __func__, __LINE__, s[i], i + 1u);
			}
			return 1;
		}
		else if (j == '-')
			continue;
		else
			s[i] = (char)(lookup_table[(unsigned int)(j) - DNA_BEGIN] + DNA_BEGIN);
	}
	return 0;
}
