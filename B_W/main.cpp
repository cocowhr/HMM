#include "viterbi.h"


//bool compstring(const seq a, const seq b) {
//	return a.first<b.first;
//}
bool compnum(const seq a, const seq b) {
	return a.second>b.second;
}
class Seq_eq
{
public:
	Seq_eq(const char& ss):s(ss){}
	bool operator() (const seq& c) const
	{
		return c.first == s;
	}
private:
	char s;
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

	return str;
}
string int2str(int n) {

	char t[24];
	int i = 0;

	while (n) {
		t[i++] = (n % 10) + '0';
		n /= 10;
	}
	t[i] = 0;
	string ret=(string)_strrev(t);
	if (ret.length()<1)
	{
		ret="0";
	}
	return ret;
}
string&   replace_all_distinct(string& str,const   string&   old_value,const   string&   new_value)     
{     
	for(string::size_type   pos(0);   pos!=string::npos;   pos+=new_value.length())   {     
		if(   (pos=str.find(old_value,pos))!=string::npos   )     
			str.replace(pos,old_value.length(),new_value);     
		else   break;     
	}     
	return   str;     
}  
double Mtest(string s,vector<seq>*sequences,double** A,int N,double* pi)
{
	int slen=s.length();
	for(int i=0;i<slen;i++)
	{
		if(s[i]>='A'&&s[i]<='Z')
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
			string temp1="";
			temp1+=s[i];
			string temp2=int2str(sid); 
			replace_all_distinct(s,temp1,temp2);
		}
	}
	int sup=0;
	for(int i=0;i<=slen-4;i++)
	{
		string ss=s.substr(i,4);
		int *v=Viterbi(A,N,ss[0]-'0',pi);
		int acc=0;
		for(int i=0; i<SEQUENCELEN; ++i)
		{
			if(v[i]==ss[i]-'0')
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

int main()
{
	string items[ITEMNUM];
	items[0]="ABCD";
	items[1]="BCDA";
	items[2]="CDAB";
	items[3]="DABC";
	items[4]="ABDE";
	items[5]="BEFG";
	items[6]="ACDE";
	items[7]="BEFC";
	items[8]="ABEF";
	items[9]="ACDC";
	map<char,int>itemmap;
	int total=0;
	for(int i=0; i<ITEMNUM; i++)
	{
		int len=items[i].length();
		for(int j=0; j<len; j++)
		{
			total++;
			char word=items[i][j];
			map<char,int>::iterator it;
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
	map<char,int>::iterator iter;
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
						items[k][m]='*';
					}
				}
			}
			minorpercent+=sequences[minor].second;
			sequences.pop_back();
		}
		sequences.push_back(make_pair('*',minorpercent));
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

	string s="GABCDA";
	double *pi=new double[num];
	for(int i=0;i<num;i++)
	{
		pi[i]=log(sequences[i].second);
	}
	cout<<Mtest(s,&sequences,A,num,pi)<<endl;

}