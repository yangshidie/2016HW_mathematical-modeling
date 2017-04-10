#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<algorithm>
#define N 12
using namespace std;



vector<string> keyrs;



void readkeyrs()
{
	FILE *fp;
	if ((fp = fopen("keyrs.dat", "r")) == NULL)
	{
		printf("Open file ERROR !\n");
		exit(0);
	}

	char *cc;
	cc = (char *)malloc(N*sizeof(char));
	fgets(cc, N, fp);
	keyrs.push_back(cc);
	while (!feof(fp))
	{
		fgets(cc, N, fp);
		keyrs.push_back(cc);
	}
	keyrs.pop_back();

	fclose(fp);
	free(cc);

	for (int i = 0; i<keyrs.size(); i++)
		printf("%s", keyrs[i].c_str());
	printf("\n");
}


bool cmprs(char *rs_c)
{
	for (int i = 0; i < keyrs.size(); i++)
	{
		if (!strcmp(keyrs[i].c_str(), rs_c))
			return true;
	}
	return false;
}


void getgefilename(int geno)
{
	char *fname;
	fname = (char *)malloc(50 * sizeof(char));
	sprintf(fname, "gene_%u.dat", geno);
	FILE *fp;
	if ((fp = fopen("keygenoname.dat", "w+")) == NULL)
	{
		printf("Open file ERROR !\n");
		exit(0);
	}

	fprintf(fp, "%s\n", fname);
	printf("%s\n", fname);
	free(fname);
	fclose(fp);
}


void main()
{
	readkeyrs();
	char *f_name;
	f_name = (char *)malloc(50 * sizeof(char));
	for (int i = 0; i<300; i++)
	{
		sprintf(f_name, "gene_info\\gene_%u.dat", i + 1);

		FILE *fp;
		if ((fp = fopen(f_name, "r")) == NULL)
		{
			printf("Open file ERROR !%d\n", i);
			exit(0);
		}
		int fte = 0;
		fseek(fp, 0, 2);
		fte = ftell(fp);
		fseek(fp, 0, 0);
		char *rs_cmp;
		rs_cmp = (char *)malloc(N*sizeof(char));
		fgets(rs_cmp, N, fp);
		while (fte != ftell(fp))
		{
			int sl = strlen(rs_cmp);
			for (int j = 0; ((rs_cmp[j]>'9') || (rs_cmp[j]<'0')) && (rs_cmp[j] != 'r') && (rs_cmp[j] != 's'); j++)
			{
				rs_cmp[j] = '\0';
			}
			if (cmprs(rs_cmp) == true)
				getgefilename(i + 1);
			fgets(rs_cmp, N, fp);
		}

		fclose(fp);
	}
	free(f_name);

}