#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXCHROM 11

char chrom[MAXCHROM + 1][8];
int nchrom;
double maxx[MAXCHROM];
double maxy[MAXCHROM];

#define MAXNXY 512

double xy[MAXCHROM][MAXNXY + 1][2];
int nxy[MAXCHROM];

int xycmp(const void *xy1, const void *xy2)
{
	double *xyd1, *xyd2;

	xyd1 = (double *) xy1;
	xyd2 = (double *) xy2;
	if (xyd1[0] == xyd2[0])
		return 0;
	if (xyd1[0] < xyd2[0])
		return -1;
	return 1;
}

double lr(int c, int x1, int x2, double *dev, double *y1, double *y2)
{
	double x, y, nd, sx, sy, sx2, sxy, denom, a, b, d, maxd;
	int xn;
	double w[MAXNXY + 1];

	nd = (double) (x2 - x1 + 1);
	for (xn = x1, sx = 0.0, sy = 0.0, sx2 = 0.0, sxy = 0.0; xn <= x2; xn++)
	{
		x = xy[c][xn][0];
		y = xy[c][xn][1];
		sx += x;
		sy += y;
		sx2 += x * x;
		sxy += x * y;
	}
	denom = nd * sx2 - sx * sx;
	a = (nd * sxy - sx * sy);
	if (fabs(denom) > a * 1e-10)
	{
		a /= denom;
		b = (sy - a * sx) / nd;
		for (xn = x1, maxd = 0.0; xn <= x2; xn++)
		{
			d = fabs(xy[c][xn][1] - (a * xy[c][xn][0] + b));
			if (d > maxd)
				maxd = d;
		}
		*dev = maxd;
		*y1 = a * xy[c][x1][0] + b;
		*y2 = a * xy[c][x2][0] + b;
		return a;
	}
	return 0.0 / 0.0;
}

void main()
{
	FILE *r;
	int ch;
	double slope, slopemax, maxslope, dev, maxdev;
	double int1, intmax, y1, y2, y1max, y2max;
	int chromn, xxn;
	int x1, x2, x1max, x2max;
	char SNP_ID[64], chr[8], chr0[8];
	double x, endx, y;

	maxslope = 2.5e-2;
	maxdev = 2.0;
	if ((r = fopen("recomb.tsv", "r")))
	{
		while ((ch = getc(r)) != EOF && ch != '\n')
			;
		for (nchrom = 0, chr0[0] = '\0';
				fscanf(r, "%s %s %lf %lf %lf",
								SNP_ID, chr, &x, &endx, &y) == 5;
				)
		{
			if (SNP_ID[0] != '#')
			{
				x /= 1e6;
				if (strcmp(chr, chr0))
				{
					if (nchrom < MAXCHROM)
					{
						strcpy(chr0, chr);
						chromn = nchrom;
						strcpy(chrom[nchrom++], chr);
						maxx[chromn] = 0.0;
						maxy[chromn] = 0.0;
					}
				}
				xy[chromn][nxy[chromn]][0] = x;
				xy[chromn][nxy[chromn]][1] = y;
				if (x > maxx[chromn])
					maxx[chromn] = x;
				if (y > maxy[chromn])
					maxy[chromn] = y;
				nxy[chromn]++;
			}
		}
		for (chromn = 0; chromn < nchrom; chromn++)
		{
			printf("nxy[%d] = %d\n", chromn, nxy[chromn]);
			qsort(xy[chromn], nxy[chromn], sizeof(double [2]), xycmp);
		}
		printf("slope limit %lg\n", maxslope);
		for (chromn = 0; chromn < nchrom; chromn++)
		{
			for (x1 = 0, intmax = 0.0; x1 < nxy[chromn] - 3; x1++)
			{
				for (x2 = x1 + 3; x2 < nxy[chromn]; x2++)
				{
					slope = lr(chromn, x1, x2, &dev, &y1, &y2);
					if (!maxslope ||
							(finite(slope) &&
								0.0 < slope && fabs(slope) <= maxslope && dev <= maxdev))
					{
						int1 = xy[chromn][x2][0] - xy[chromn][x1][0];
						if (!intmax || int1 > intmax)
						{
							x1max = x1;
							x2max = x2;
							y1max = y1;
							y2max = y2;
							intmax = int1;
							slopemax = slope;
						}
					}
				}
			}
			printf("%s: x1 = %lg, y1 = %lg, x2 = %lg, y2 = %lg\n",
											chrom[chromn], xy[chromn][x1max][0], y1max,
											xy[chromn][x2max][0], y2max);
		}
	}
}
