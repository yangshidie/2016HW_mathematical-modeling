#include<stdio.h>
#include<stdio.h>
#include<algorithm>
#include<vector>
#define NP 25
#define TR 424
#define N 100000
#define ER 10
#define NU 99
#define DVALUE 0.12
using namespace std;



vector<int> hqu;
vector<int> uhqu;
vector<int> rs_qu;
vector<float>rs_pe;
short h_rs[TR][1000];        //无病位点信息
short uh_rs[TR][1000];			//有病位点信息
bool orph[1000][10];    //原始数据数组



void getexchqu()
{

	FILE *fp;
	if ((fp = fopen("multi_phenos.txt", "r")) == NULL)
	{
		printf("Open file ERROR !\n");
		exit(0);
	}


	char *cc;
	cc = (char *)malloc(NP*sizeof(char));
	fseek(fp, 0, 2);
	int fte = ftell(fp);
	fseek(fp, 0, 0);
	int pe_no = 0;
	char c;
	fgets(cc, NP, fp);
	while (ftell(fp) != fte)
	{
		for (int i = 0; i < 10; i++)
		{
			c = cc[i * 2];
			switch (c)
			{
			case'0':
				orph[pe_no][i] = 0;
				break;
			case'1':
				orph[pe_no][i] = 1;
				break;
			default:
				break;
			}
		}
		pe_no++;
		fgets(cc, NP, fp);
	}

	free(cc);
	fclose(fp);
}



void test1()
{
	int no_1 = 0;
	int no_0 = 0;
	for (int i = 0; i < 1000; i++)
	{
		bool bo_1 = 1;
		bool bo_0 = 0;
		for (int j = 0; j < 10; j++)
		{
			bo_1 = bo_1 & orph[i][j];
			bo_0 = bo_0 | orph[i][j];
		}
		if (bo_1 == 1)
		{
			printf("%d : all 1\n", i);
			no_1++;
		}
		if (bo_0 == 0)
		{
			printf("%d : all 0\n", i);
			no_0++;
		}
	}
	printf("ALL 0 :%d\n", no_1);
	printf("ALL 1 :%d\n", no_0);

}


void test2()
{
	FILE *fp;
	if ((fp = fopen("geno_count.dat", "w")) == NULL)
	{
		printf("Open file ERROR !\n");
		exit(0);
	}


	int a[11];
	for (int i = 0; i < 11; i++)
		a[i] = 0;
	for (int i = 0; i < 1000; i++)
	{
		int no = 0;
		for (int j = 0; j < 10; j++)
		{
			if (orph[i][j] == true)
				no++;
		}
		a[no]++;
	}
	for (int i = 0; i < 11; i++)
	{
		printf("%d :%d\n", i, a[i]);
		fprintf(fp, "%d  %d\n", i, a[i]);
	}


	fclose(fp);
}


void getkeypeno()                 //已验证正确
{
	for (int i = 0; i < 1000; i++)
	{
		int no = 0;
		for (int j = 0; j < 10; j++)
		{
			if (orph[i][j] == true)
				no++;
		}
		if ((no == 0) || (no == 1) || (no == 2))
		{
			hqu.push_back(i);
		}
		if ((no == 10) || (no == 9) || (no == 8))
		{
			uhqu.push_back(i);
		}
	}
}


void readrs(int cx)       //已验证正确
{
	FILE *fp;
	if ((fp = fopen("genotypebi.dat", "r")) == NULL)
	{
		printf("Open file ERROR !/n");
		exit(0);
	}


	int i_nu = 0, pos;
	char *cc, tempc;
	cc = (char*)malloc(N*sizeof(char));
	fgets(cc, N, fp);
	while ((!feof(fp)) && (i_nu < 1000))
	{
		fgets(cc, N, fp);

		if (find(hqu.begin(), hqu.end(), i_nu) != hqu.end())      //该行在全健康序列中
		{
			pos = 0;
			while (i_nu != hqu[pos])        //求该人在健康序列位置
				pos++;
			if (pos >= TR)
				continue;
			for (int i = 0; i < 1000; i++)
			{
				tempc = cc[i * 2 + cx * 2000];
				if ((tempc >= '0') && (tempc <= '9'))
					h_rs[pos][i] = tempc - '0';
				else
					h_rs[pos][i] = ER;
			}
			if (cx >= 9)
			{
				for (int j = 445; j < 1000; j++)
					h_rs[pos][j] = NU;
			}

		}


		if (find(uhqu.begin(), uhqu.end(), i_nu) != uhqu.end())      //该行在全有病序列中
		{
			pos = 0;
			while (i_nu != uhqu[pos])        //求该人在全有病序列位置
				pos++;
			if (pos >= TR)
				continue;
			for (int i = 0; i < 1000; i++)
			{
				tempc = cc[i * 2 + cx * 2000];
				if ((tempc >= '0') && (tempc <= '9'))
					uh_rs[pos][i] = tempc - '0';
				else
					uh_rs[pos][i] = ER;
			}
			if (cx >= 9)
			{
				for (int j = 445; j < 1000; j++)
					uh_rs[pos][j] = NU;
			}
		}
		i_nu++;
	}

	free(cc);
	fclose(fp);
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
	for (i = 0; i < TR; i++)
		swh(uh_c, uh_rs[i][rs_no % 1000]);

	for (i = 0; i<TR; i++)
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


void rers()
{
	int cx = 0;
	for (cx = 0; cx < 10; cx++)
	{
		readrs(cx);
		if (cx < 9)
		{
			for (int rs_no = 0; rs_no < 1000; rs_no++)
			{
				rs_count(rs_no + cx * 1000, cx);
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
	printf("Sample size: N=%d\nDVALUE =%.4f\n\n", TR ,DVALUE);
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
				printf("%s\t\t%f\n", temps,  rs_pe[veno]);
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


void main()
{
	getexchqu();         //读取multi_pheno
	getkeypeno();        //求关键人群
	rers();

	outrst();
}

