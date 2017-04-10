#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<algorithm>
#include<vector>
#define N 100000
#define ER 10
#define DVALUE 0.115
#define TR 500
#define NU 99
using namespace std;



short uh_rs[500][1000];
short h_rs[500][1000];
vector<int> rs_qu;
vector<float>rs_pe;



void readrs(int cx)
{
	FILE *fp;
	if ((fp = fopen("genotypebi.dat", "r")) == NULL)
	{
		printf("Open file ERROR !\n");
		exit(0);
	}

	int i_nu = 0;
	char *cc, tempc;
	cc = (char*)malloc(N*sizeof(char));
	fgets(cc, N, fp);
	while ((!feof(fp)) && (i_nu<1000))
	{
		fgets(cc, N, fp);

		if (i_nu<500)
		{
			for (int i = 0; i<1000; i++)
			{

				tempc = cc[i * 2 + cx * 2000];
				if ((tempc >= '0') && (tempc <= '9'))
					h_rs[i_nu][i] = tempc - '0';
				else
					h_rs[i_nu][i] = ER;
			}
			if (cx >= 9)
			{
				for (int j = 445; j<1000; j++)
					h_rs[i_nu][j] = NU;
			}
		}
		else
		{
			for (int i = 0; i<1000; i++)
			{
				tempc = cc[i * 2 + cx * 2000];
				if ((tempc >= '0') && (tempc <= '9'))
					uh_rs[i_nu - 500][i] = tempc - '0';
				else
					uh_rs[i_nu - 500][i] = ER;
			}
			if (cx >= 9)
			{
				for (int j = 445; j<1000; j++)
					uh_rs[i_nu - 500][j] = NU;
			}
		}


		i_nu++;
	}


	free(cc);
	fclose(fp);
}


bool getrs_no(float pe, int rs_no)
{
	if (pe >= DVALUE)
	{
		rs_qu.push_back(rs_no);
		rs_pe.push_back(pe);
		return true;
	}
	return false;
}


void swh(int *u, int uh_rs)
{
	switch (uh_rs)
	{
	case 0:
		u[0]++;
		break;
	case 1:
		u[1]++;
		break;
	case 2:
		u[2]++;
		break;
	case 3:
		u[3]++;
		break;
	case 4:
		u[4]++;
		break;
	case 5:
		u[5]++;
		break;
	case 6:
		u[6]++;
		break;
	case 7:
		u[7]++;
		break;
	case 8:
		u[8]++;
		break;
	case 9:
		u[9]++;
		break;
	case 10:
		u[10]++;
		break;
	default:
		break;
	}
}


void rs_count(int rs_no, int cx)
{
	int h_c[11], h_sum;
	int uh_c[11], uh_sum;
	int i;

	float h_pe[11];
	float uh_pe[11];


	for (i = 0; i<11; i++)
	{
		h_c[i] = 0;
		uh_c[i] = 0;
		h_pe[i] = 0;
		uh_pe[i] = 0;
	}
	for (i = 0; i < TR; i++)          //TR 个建立模型，剩下的检验
		swh(uh_c, uh_rs[i][rs_no % 1000]);

	for (i = 0; i<TR; i++)          //TR建立模型，剩下100检验
		swh(h_c, h_rs[i][rs_no % 1000]);

	h_sum = 0;
	uh_sum = 0;
	for (i = 0; i<11; i++)
	{
		h_sum += h_c[i];
		uh_sum += uh_c[i];
	}
	if (h_sum != TR)
		printf("%d SUM ERROR\n", h_sum);
	if (uh_sum != TR)
		printf("%d	UH SUM ERROR\n", uh_sum);

	for (i = 0; i < 11; i++)
	{
		h_pe[i] = (float)h_c[i] / (float)TR;
		uh_pe[i] = (float)uh_c[i] / (float)TR;
	}


	float pemax = abs(h_pe[0] - uh_pe[0]);
	for (i = 1; i < 11; i++)
	{
		if (pemax<abs(h_pe[i] - uh_pe[i]))
			pemax = abs(h_pe[i] - uh_pe[i]);
	}
	bool flag = getrs_no(pemax, rs_no);

}


int findno(int nu)
{
	for (int i = 0; i<rs_qu.size(); i++)
	{
		if (nu == rs_qu[i])
			return i;
	}
	return N;
}


void outrst()
{
	FILE *fp;
	if ((fp = fopen("genotypebi.dat", "r")) == NULL)
	{
		printf("Open file ERROR !/n");
		exit(0);
	}


	FILE *out;
	if ((out = fopen("keyrs_pe.dat", "w")) == NULL)
	{
		printf("Open file ERROR !/n");
		exit(0);
	}


	FILE *outrs;
	if ((outrs = fopen("keyrs.dat", "w")) == NULL)
	{
		printf("Open file ERROR !/n");
		exit(0);
	}


	char ch[16] = { "Sample size: N=" };
	fputs(ch, out);
	fprintf(out, "%d\n", TR);
	char ch1[16] = { "DVALUE=" };
	fputs(ch1, out);
	fprintf(out, "%.4f\n", DVALUE);
	printf("Sample size: N=%d\nDVALUE=%.4f\n\n", TR, DVALUE);
	char *cc, temps[11];
	cc = (char*)malloc(N*sizeof(char));
	fgets(cc, N, fp);

	int nu = 0, t_i = 0;

	for (int i = 0; cc[i] != NULL; i++)
	{

		if (cc[i] == 'r')
		{
			temps[0] = cc[i];
			t_i = 0;
		}
		else if ((cc[i] != 's') && ((cc[i]<'0') || (cc[i]>'9')))
		{
			t_i++;
			temps[t_i] = '\0';
			int veno = findno(nu);
			if (veno != N)
			{
				fprintf(out, "%s\t%f\n", temps, rs_pe[veno]);
				fprintf(outrs, "%s\n", temps);
				printf("%s\t\t%f\n", temps, rs_pe[veno]);
			}
			nu++;
		}
		else
		{
			t_i++;
			temps[t_i] = cc[i];
		}
	}


	free(cc);
	fclose(fp);
	fclose(out);
	fclose(outrs);
}


int main()
{
	int cx = 0;    //每次1000个循环9次

	for (cx = 0; cx<10; cx++)
	{
		readrs(cx);
		if (cx<9)
		{
			for (int rs_no = 0; rs_no<1000; rs_no++)       //rs_no位点
			{
				rs_count(rs_no + cx * 1000, cx);//统计
			}
		}
		else
		{
			for (int rs_no = 0; rs_no<445; rs_no++)       //rs_no位点
			{
				rs_count(rs_no + cx * 1000, cx);//基因统计
			}
		}
	}

	outrst();
	return 0;
}




