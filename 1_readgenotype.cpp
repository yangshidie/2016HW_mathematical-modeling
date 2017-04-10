#include<stdio.h>
#include<stdlib.h>
#include<string>
#define N 100000
using namespace std;



char swtch(char c1, char c2)
{
	char ctem = 'X';
	switch (c1)
	{
	case 'A':
		switch (c2)
		{
		case 'A':
			ctem = '0';
			break;
		case 'C':
			ctem = '1';
			break;
		case 'G':
			ctem = '2';
			break;
		case 'T':
			ctem = '3';
			break;
		default:
			break;
		}
		break;

	case 'C':
		switch (c2)
		{
		case 'A':
			ctem = '1';
			break;
		case 'C':
			ctem = '4';
			break;
		case 'G':
			ctem = '5';
			break;
		case 'T':
			ctem = '6';
			break;
		default:
			break;
		}
		break;

	case 'G':
		switch (c2)
		{
		case 'A':
			ctem = '2';
			break;
		case 'C':
			ctem = '5';
			break;
		case 'G':
			ctem = '7';
			break;
		case 'T':
			ctem = '8';
			break;
		default:
			break;
		}
		break;

	case 'T':
		switch (c2)
		{
		case 'A':
			ctem = '3';
			break;
		case 'C':
			ctem = '6';
			break;
		case 'G':
			ctem = '8';
			break;
		case 'T':
			ctem = '9';
			break;
		default:
			break;
		}
		break;
	default:
		break;
	}
	return ctem;
}


int main()
{
	FILE *in, *out;
	char infile[14], outfile[16];

	strcpy(infile, "genotype.dat");
	if ((in = fopen("genotype.dat", "r")) == NULL)
	{
		printf("Open file ERROR !\n");
		exit(0);
	}

	strcpy(outfile, "genotypebi.dat");
	if ((out = fopen("genotypebi.dat", "w")) == NULL)
	{
		printf("Open file ERROR !\n");
		exit(0);
	}

	int ba = 0, i = 0;
	string s;
	char *cc, c, te;
	cc = (char*)malloc(N*sizeof(char));

	fgets(cc, N, in);
	fputs(cc, out);
	int col = 0;
	while ((!feof(in)) && (col<1000))
	{
		fgets(cc, N, in);
		for (i = 0; cc[ba + i] != '\0'; i = i + 3)
		{
			te = swtch(cc[ba + i], cc[ba + i + 1]);
			fputc(te, out);
			fputc(9, out);
		}
		fputc(10, out);
		col++;
	}

	printf("Success!\n");
	free(cc);
	fclose(in);
	fclose(out);
	return 0;
}


/*

A	C	G	T

A	0	1	2	3

C	1	4	5	6

G	2	5	7	8

T	3	6	8	9

Òì³£ X

*/