#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>

#define DEBUG 0

#define MAXGENE 8192
#define CONTLEN 63

typedef struct
{
	char con1[CONTLEN + 1];
	int bp1a, bp1b;
	char con2[CONTLEN + 1];
	int bp2a, bp2b;
} gene;

gene contig[MAXGENE];

int getgene(FILE *homol, gene *g)
{
	int ch, c;

	for (c = 0; (ch = getc(homol)) != EOF && ch != '\n' && ch != '\t'; )
	{
		if (c < CONTLEN)
			g->con1[c++] = ch;
		else
		{
			g->con1[c] = '\0';
			fprintf(stderr, "%s too long\n", g->con1);
			exit(1);
		}
	}
	g->con1[c] = '\0';
	if (ch == EOF)
		return 0;
	if (ch != '\t')
	{
		fprintf(stderr, "%s not tab terminated\n", g->con1);
		exit(1);
	}
	for (g->bp1a = 0; (ch = getc(homol)) >= '0' && ch <= '9'; )
		g->bp1a = 10 * g->bp1a + ch - '0';
	if (ch != '\t')
	{
		fprintf(stderr, "%s: bp1a not tab terminated\n", g->con1);
		exit(1);
	}
	for (g->bp1b = 0; (ch = getc(homol)) >= '0' && ch <= '9'; )
		g->bp1b = 10 * g->bp1b + ch - '0';
	if (ch != '\t')
	{
		fprintf(stderr, "%s: bp1b not tab terminated\n", g->con1);
		exit(1);
	}
	for (c = 0; (ch = getc(homol)) != EOF && ch != '\n' && ch != '\t'; )
	{
		if (c < CONTLEN)
			g->con2[c++] = ch;
		else
		{
			g->con2[c] = '\0';
			fprintf(stderr, "%s too long\n", g->con2);
			exit(1);
		}
	}
	g->con2[c] = '\0';
	if (ch != '\t')
	{
		fprintf(stderr, "%s not tab terminated\n", g->con2);
		exit(1);
	}
	for (g->bp2a = 0; (ch = getc(homol)) >= '0' && ch <= '9'; )
		g->bp2a = 10 * g->bp2a + ch - '0';
	if (ch != '\t')
	{
		fprintf(stderr, "%s: bp2a not tab terminated\n", g->con1);
		exit(1);
	}
	for (g->bp2b = 0; (ch = getc(homol)) >= '0' && ch <= '9'; )
		g->bp2b = 10 * g->bp2b + ch - '0';
	if (ch != '\n' && ch != '\t')
	{
		fprintf(stderr, "%s: bp1b not tab or newline terminated\n", g->con1);
		exit(1);
	}
	while (ch != '\n' && ch != EOF)
		ch = getc(homol);
	if (DEBUG >= 4)
	{
		printf("%s %d %d %s %d %d\n",
				g->con1, g->bp1a, g->bp1b,
				g->con2, g->bp2a, g->bp2b);
	}
	return 1;
}

int getcont(FILE *homol, int *ngene)
{
	int gotgene;

	if (DEBUG >= 4)
		printf("getcont: *ngene = %d\n", *ngene);
	if (*ngene)
	{
		strcpy(contig[0].con1, contig[*ngene].con1);
		contig[0].bp1a = contig[*ngene].bp1a;
		contig[0].bp1b = contig[*ngene].bp1b;
		strcpy(contig[0].con2, contig[*ngene].con2);
		contig[0].bp2a = contig[*ngene].bp2a;
		contig[0].bp2b = contig[*ngene].bp2b;
		*ngene = 1;
	}
	for (gotgene = 1; gotgene && *ngene < MAXGENE; )
	{
		if (DEBUG >= 4)
			printf("getcont: while: *ngene = %d\n", *ngene);
		if (*ngene < MAXGENE - 1)
		{
	 		gotgene = getgene(homol, contig + *ngene);
	 		if (!gotgene ||
							(*ngene >= 1 &&
							(strcmp(contig[*ngene].con1, contig[*ngene - 1].con1) ||
						 	strcmp(contig[*ngene].con2, contig[*ngene - 1].con2))))
				return gotgene;
			++*ngene;
		}
		else
		{
			fprintf(stderr, "Too many genes after %s\n", contig[*ngene].con1);
			exit(1);
		}
	}
	return 0;
}

void countsyn(int ngene, int *totgene, int *nsyn)
{
	int g, gs;
	int syn, csyn, csynr;

	if (DEBUG >= 2)
	{
		for (g = 0; g < ngene; g++)
		{
			printf("%s %d %d %s %d %d\n",
					contig[g].con1, contig[g].bp1a, contig[g].bp1b,
					contig[g].con2, contig[g].bp2a, contig[g].bp2b);
		}
		putchar('\n');
	}
	for (g = 1; g < ngene; g++)
	{
		if (contig[g].bp1a < contig[g - 1].bp1a)
			printf("%s out of order at %d\n", contig[0].con1, contig[g].bp1a);
	}
	for (g = 0; g < ngene; g++)
	{
		if (contig[g].bp1a > contig[g].bp1b)
		{
			if (DEBUG >= 1)
				printf("%d > %d\n", contig[g].bp1a, contig[g].bp1b);
			contig[g].bp1a = contig[g].bp1b;
		}
		if (contig[g].bp2a > contig[g].bp2b)
		{
			if (DEBUG >= 1)
				printf("%d > %d\n", contig[g].bp2a, contig[g].bp2b);
			contig[g].bp2a = contig[g].bp2b;
		}
	}
	for (g = 0, csyn = 0; g < ngene; )
	{
		syn = 1; /* start a syn */
		for (gs = g + 1; gs < ngene && contig[gs].bp2a > contig[gs - 1].bp2a; gs++)
			syn++;
		if (syn >= 3)
			csyn += syn;
		g = gs;
	}
	for (g = 0, csynr = 0; g < ngene; )
	{
		syn = 1; /* start a syn */
		for (gs = g + 1; gs < ngene && contig[gs].bp2a < contig[gs - 1].bp2a; gs++)
			syn++;
		if (syn >= 3)
			csynr += syn;
		g = gs;
	}
	if (DEBUG >= 1)
		printf("%s - %s:\n%d sytenies and %d reverse syntenies out of %d genes\n",
										contig[0].con1, contig[0].con2, csyn, csynr, ngene);
	if (DEBUG >= 2)
	{
		for (g = 0; g < ngene; g++)
		{
			printf("%s %d %d %s %d %d\n",
					contig[g].con1, contig[g].bp1a, contig[g].bp1b,
					contig[g].con2, contig[g].bp2a, contig[g].bp2b);
		}
		putchar('\n');
	}
	if (csynr > csyn)
		*nsyn += csynr;
	else
		*nsyn += csyn;
	*totgene += ngene;
}

void main(int argc, char **argv)
{
	FILE *homol;
	int ngene, gotcont, nsyn, totgene;

	puts("Filename\t# Syntenies\t# Homologs\tFraction Syntenies");
	while (--argc > 0)
	{
		if ((homol = fopen(*++argv, "r")))
		{
			for (ngene = 0, nsyn = 0, totgene = 0, gotcont = 1; gotcont; )
			{
				gotcont = getcont(homol, &ngene);
				countsyn(ngene, &totgene, &nsyn);
			}
			printf("%s\t%d\t%d\t%lf\n", *argv, nsyn, totgene, ((double) nsyn) / ((double) totgene));
			fclose(homol);
		}
		else
			fprintf(stderr, "Can't read %s\n", *argv);
	}
}
