#include<stdio.h>
#include<string.h>
const int maxn=5005;
struct Edge
{
	int s,t;
	int next;
	int vis;
}edge[1000005];
int head[maxn];
int E=0;
void add_edge(int s,int t)
{
	edge[E].s=s;
	edge[E].t=t;
	edge[E].vis=0;
	edge[E].next=head[s];
	head[s]=E++;
}
int Btype,Time,N,M;
int dfn[maxn],low[maxn],Belong[maxn];
int  st[maxn],Top;
inline int min(int a,int b){return a<b?a:b;}
void dfs(int s)
{
	int i,t,id;
	st[++Top]=s;
	dfn[s]=low[s]=++Time;
	for (i=head[s];i!=-1;i=edge[i].next)
	{	
		if(edge[i].vis)continue;	
		edge[i].vis=edge[i^1].vis=1;
		t=edge[i].t;
		if (!dfn[t])
		{
			dfs(t);
			low[s]=min(low[s],low[t]);
		}
		else low[s]=min(low[s],dfn[t]);
	}
	if(dfn[s]==low[s])
	{
		Btype++;
		do{
			t=st[Top--];
			Belong[t]=Btype;
		}while(t!=s);
	}
}
void SCC(int n)
{
	int i;
	Time=0;Btype=0;Top=0;
	memset(dfn,0,sizeof(int)*(n+1));
	for(i=1;i<=n;i++)if(!dfn[i])
		dfs(i);
}
