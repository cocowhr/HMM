#include "viterbi.h"
/********/
/*适用于较长序列的挖掘*/
/********/

//bool compstring(const seq a, const seq b) {
//	return a.first<b.first;
//}
bool compnum(const seq a, const seq b) {
	return a.second>b.second;
}
class Seq_eq
{
public:
	Seq_eq(const int& ss):s(ss){}
	bool operator() (const seq& c) const
	{
		return c.first == s;
	}
private:
	int s;
};
///*注意：当训练序列为ABCD时
//探查以B开始长度为4的序列会出现错误的后果
//TODO:fix
int* Viterbi(double **A,int N, int s,double* pi)
{ 
	int T=SEQUENCELEN;
	int* str=new int[T];
	str[0]=s;
	double **delta=new double*[T];
	for(int i=0;i<T;i++)
	{
		delta[i]=new double[N]();
	}
	int **psi=new int*[T];
	for(int i=0;i<T;i++)
	{
		psi[i]=new int[N]();
	}
	int 	i, j;	/* 状态下标 */
	int  	t;	/* 时间下标 */

	int	minvalind;
	double	minval, val;
	/* 1. 初始化  */
	for (i = 0; i < N; i++)
	{
		delta[1][i] = pi[i] + A[s][i];
		psi[1][i] = s;
	}
	/* 2. 递归 */
	for (t = 2; t < T; t++)
	{
		for (j = 0; j < N; j++)
		{
			minval = delta[t-1][0]+(A[0][j]);
			minvalind = 0;
			for (i = 1; i < N; i++)
			{
				val = delta[t-1][i]+(A[i][j]);
				if (val > minval)
				{
					minval = val;
					minvalind = i;
				}
			}
			delta[t][j] = minval;
			psi[t][j] = minvalind;
		}
	}

	/* 3. 终止 */
	double pprob = delta[T-1][0];
	str[T-1]=0;
	for (i = 1; i <N; i++)
	{
		if (delta[T-1][i] > pprob)
		{
			pprob = delta[T-1][i];
			str[T-1] = i;
		}
	}
	/* 4. Path (state sequence) backtracking */
	for (t = T - 2; t >= 1; t--)
		str[t] = psi[t+1][str[t+1]];
	for(int i=0;i<T;i++)
	{
		cout<<str[i];
	}
	cout<<endl;
	return str;
}
double Mtest(int s[],vector<seq>*sequences,double** A,int N,double* pi)
{
	int slen=4;
	for(int i=0;i<slen;i++)
	{
		int sid=0;
		vector<seq>::iterator result = find_if( sequences->begin( ), sequences->end( ),Seq_eq(s[i])); 
		if(result != sequences->end())
		{
			sid=result-sequences->begin();
		}
		else
		{
			sid=sequences->end()-sequences->begin()-1;
		}
		s[i]=sid;
	}
	int sup=0;
	for(int i=0;i<=slen-4;i++)
	{
		int ss[4];
		ss[0]=s[i];
		ss[1]=s[i+1];
		ss[2]=s[i+2];
		ss[3]=s[i+3];
		int *v=Viterbi(A,N,ss[0],pi);
		int acc=0;
		for(int i=0; i<SEQUENCELEN; ++i)
		{
			/*cout<<ss[i]<<endl;*/
			if(v[i]==ss[i])
			{
				acc++;
			}
		}
		if(((double)(acc)/SEQUENCELEN)>TV)
		{
			sup++;
		}
	}
	return (double)sup/(slen-SEQUENCELEN+1);
}
double Mtestacc(int s[],vector<seq>*sequences,double** A,int N,double* pi)
{
	for(int i=0;i<4;i++)
	{
		int sid=0;
		vector<seq>::iterator result = find_if( sequences->begin( ), sequences->end( ),Seq_eq(s[i])); 
		if(result != sequences->end())
		{
			sid=result-sequences->begin();
		}
		else
		{
			sid=sequences->end()-sequences->begin()-1;
		}
		s[i]=sid;
	}
	int sup=0;
	int *v=Viterbi(A,N,s[0],pi);
	int acc=0;
	for(int i=0; i<SEQUENCELEN; ++i)
	{
		/*cout<<ss[i]<<endl;*/
		if(v[i]==s[i])
		{
			acc++;
		}
	}
	return (double)(acc)/SEQUENCELEN;
}
int main()
{
	int items[ITEMNUM][ITEMLEN];
	mysql_init(&mysql);  
	mysql_real_connect(&mysql, "localhost", "root", "root", "genet", 3306, NULL, 0); 
	strSQL.Format("SELECT * FROM seq2");
	mysql_query(&mysql, strSQL);
	result = mysql_store_result(&mysql);  
	int n=0;
	while ((row = mysql_fetch_row(result))!=NULL)  
	{  
		for(int i=0;i<ITEMLEN;i++)
		{
			items[n][i]=atoi(row[i+1]);
		}
		n++;
	}  
	//items[0][0]=11;
	//items[0][1]=12;
	//items[0][2]=13;
	//items[0][3]=14;//abcd
	//items[1][0]=12;
	//items[1][1]=13;
	//items[1][2]=14;
	//items[1][3]=11;//bcda
	//items[2][0]=13;
	//items[2][1]=14;
	//items[2][2]=11;
	//items[2][3]=12;//cdab
	//items[3][0]=14;
	//items[3][1]=11;
	//items[3][2]=12;
	//items[3][3]=13;//dabc
	//items[9][0]=11;
	//items[9][1]=12;
	//items[9][2]=14;
	//items[9][3]=15;//abde
	//items[4][0]=12;
	//items[4][1]=15;
	//items[4][2]=16;
	//items[4][3]=17;//befg
	//items[5][0]=11;
	//items[5][1]=13;
	//items[5][2]=14;
	//items[5][3]=15;//acde
	//items[6][0]=12;
	//items[6][1]=15;
	//items[6][2]=16;
	//items[6][3]=13;//befc
	//items[7][0]=11;
	//items[7][1]=12;
	//items[7][2]=15;
	//items[7][3]=16;//abef
	//items[8][0]=11;
	//items[8][1]=13;
	//items[8][2]=14;
	//items[8][3]=13;//acdc
	/*items[0]="ABCD";
	items[1]="BCDA";
	items[2]="CDAB";
	items[3]="DABC";
	items[4]="ABDE";
	items[5]="BEFG";
	items[6]="ACDE";
	items[7]="BEFC";
	items[8]="ABEF";
	items[9]="ACDC";*/
	map<int,int>itemmap;
	int total=0;
	for(int i=0; i<ITEMNUM; i++)
	{
		for(int j=0; j<ITEMLEN; j++)
		{
			total++;
			int word=items[i][j];
			map<int,int>::iterator it;
			it=itemmap.find(word);
			if(it==itemmap.end())
			{
				itemmap[word]=1;
			}
			else
			{
				itemmap[word]++;
			}
		}
	}
	vector<seq> sequences;
	map<int,int>::iterator iter;
	for(iter=itemmap.begin();iter!=itemmap.end();iter++)
	{
		sequences.push_back(make_pair(iter->first,(double)iter->second/total));
	}
	itemmap.clear();
	sort(sequences.begin(),sequences.end(),compnum);
	double totalprecent=0;
	int num=sequences.size();
	int major=0;
	for(;totalprecent<=0.9&&major<num;major++)//TODO:totalprecent应为0.99 为了测试现在为0.9 小于0.01的都为少数状态
	{
		totalprecent+=sequences[major].second;
	}
	if(major<num)
	{
		double minorpercent=0;
		for(int minor=num-1;minor>=major;minor--)
		{
			for(int k=0;k<ITEMNUM;k++)
			{
				for(int m=0;m<ITEMLEN;m++)
				{
					if(sequences[minor].first==items[k][m])
					{
						items[k][m]=-1;
					}
				}
			}
			minorpercent+=sequences[minor].second;
			sequences.pop_back();
		}
		sequences.push_back(make_pair(0,minorpercent));
	}
	num=sequences.size();
	double** A=new double* [num];
	double* B=new double[num]();
	for(int i=0;i<num;i++)
	{
		for(int k=0;k<ITEMNUM;k++)
		{
			for(int m=0;m<ITEMLEN-1;m++)
			{
				if(sequences[i].first==items[k][m])
				{
					B[i]++;
				}
			}
		}
	}
	for(int i=0;i<num;i++)
	{
		A[i]=new double[num](); 
	}
	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			for(int k=0;k<ITEMNUM;k++)
			{
				for(int m=0;m<ITEMLEN-1;m++)
				{
					if(sequences[i].first==items[k][m])
					{
						if(sequences[j].first==items[k][m+1])
						{
							A[i][j]++;
						}
					}
				}
			}
		}
	}
	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			if(B[i]!=0)
			{
				A[i][j]=log(A[i][j]/B[i]);
			}
			else
			{
				A[i][j]=0;
				A[i][j]=log(A[i][j]);
			}
		}
	}

	//int s[4]={4,5,9,1};//假设命令为单数 之后会用数组代替
	double *pi=new double[num];
	for(int i=0;i<num;i++)
	{
		pi[i]=log(sequences[i].second);
	}
	int *v=Viterbi(A,num,0,pi);
	for(int i=0;i<4;i++)
		cout<<sequences[v[i]].first<<endl;
	//cout<<Mtestacc(s,&sequences,A,num,pi)<<endl;

}