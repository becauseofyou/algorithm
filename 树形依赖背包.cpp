/*
The most powerful force
hdu 3593
*/
#include<cstdio>
#include<cstring>
#include<set>
#include<string>
#include<iostream>
#include<cmath>
#include<vector>
#include<map>
#include<stack>
#include<climits>
#include<numeric>
#include<time.h>
#include<queue>
#include<cstdlib>
#include<algorithm>
using namespace std;
#define lowbit(x) ((x)&(-(x)))
#define sqr(x) ((x)*(x))
#define PB push_back
#define MP make_pair
#define foreach(e,x) for(vector<int>::iterator e=x.begin();e!=x.end();++e)
#define For(i,a,b) for(int i=a;i<=b;i++)
typedef long long lld;
typedef vector<int> VI;
typedef vector<string> VS;
typedef pair<int,int> PII;
const int maxn  =  200010;
const int inf = ~0u>>2;
inline void Max(int &a,int b){	if(b>a) a=b;}
inline void Min(int &a,int b){  if(b<a) a=b;}
int head[maxn], pnt[maxn] , nxt[maxn] ,E;
int val[maxn] , cost[maxn] ;
bool leaf[maxn];
void add(int a,int b)
{
	pnt[E] = b;
	nxt[E] = head[a];
	head[a] = E++;
}
int dp[510][10010];
void dfs(int u,int cap)
{
	for(int i=head[u];i!=-1;i=nxt[i])
	{
		int now = pnt[i];
		if(!leaf[now]) {
			for(int j=cap-cost[now];j>=0;j--)	dp[now][j] = dp[u][j] + val[now];
			dfs(now,cap-cost[now]);
			for(int j=cap-cost[now];j>=0;j--)   
			{
				Max(dp[u][j+cost[now]],dp[now][j]);
			}
		} else {
			for(int j=cap-cost[now];j>=0;j--)
			{
				Max(dp[u][j+cost[now]],dp[u][j] + val[now]);
			}
		}
	}
}
int main()
{
	int n,m;
	while(~scanf("%d%d",&n,&m))
	{
		int p;
		E=0; val[0] = cost[0] = 0;
		fill(leaf,leaf+n+1,true);
		fill(head,head+n+1,-1);
		For(i,1,n)
		{
			scanf("%d%d%d",&cost[i],&val[i],&p);
			if(p == i) p = 0;
			add(p,i);
			leaf[p]=false;
		}
		fill(dp[0],dp[0]+m+1,0);
		dfs(0,m);
		printf("%d\n",dp[0][m]);
	}
	return 0;
}
