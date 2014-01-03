//hdu 3686
#include <cstdio>
#include <algorithm>
#include <vector>
#include <cstring>
using std::vector;
const int N=50000;
const int M=1000000;
vector<int> edge[N];
struct Edge{
  int id;
  int s,t,next;
}list[M];
int head[N],E,dfn[N],low[N], stk[M], tot, n,nn, ok[N], belong[M],child, Btype,tdfn,ve[M];
void init(){
  memset(ve,false,sizeof(ve));
  memset(ok,false,sizeof(ok));
  tdfn = 0;
  Btype = 0;
  memset(dfn,0,sizeof(dfn));
}
void add(int u,int v,int _id){
  list[E].id = _id;
  list[E].t = v;
  list[E].s = u;
  list[E].next = head[u];
  head[u] = E++;
}
void dfs(int u,int f){
  low[u] = dfn[u] = ++tdfn;
  for(int i = head[u]; i != -1; i = list[i].next){
    if(ve[list[i].id])  continue;
    int v = list[i].t;
    ve[list[i].id] = true;
    stk[++tot] = i;
    if(!dfn[v]) {    
      if(u == f) child ++;//f is root
      dfs(v,f);
      if(low[v] < low[u]) low[u] = low[v];
      if(low[v] >= dfn[u]){
        if(u!=f) ok[u] = true; // u is gedian 
        ++Btype;
        int id;
        do {
          id = stk[tot--];
          belong[id/2+1] = Btype;
        }while(list[id].s!=u);
      }
    }
    else if(dfn[v] < low[u]) low[u] = dfn[v];
  }
}
void SCC(){
  init();
  for(int i = 1; i <= n; i++) {
    child = 0; tot = 0;
    if(!dfn[i]) {
      dfs(i,i);
      ok[i] = false;
      if(child > 1) ok[i] = true;
    }
  }
  nn = 0;
  for(int i = 1; i <= n; i++) if(ok[i]) nn++;
  nn += Btype;
  for(int i = 0; i < N; i++) edge[i].clear();
  for(int i = 1,now=0; i <= n; i++)if(ok[i]){
    ++now;
    for(int j = head[i]; j != -1; j = list[j].next){
      edge[now+Btype].push_back(belong[j/2+1]);
      edge[belong[j/2+1]].push_back(now+Btype);
    }
  }
}
int f[18][N];
bool vis[N];
int dep[N];
void DFS(int u,int fa) {
  vis[u] = true; dep[u] = dep[fa] + 1;
  for(int i = 0; i < edge[u].size(); i++){
    int v = edge[u][i];
    if(vis[v]) continue;
    f[0][v] = u;
    for(int j = 1; j < 18; j++) f[j][v] = f[j-1][f[j-1][v]];
    DFS(v,u);
  }
}
int LCA(int a,int b) {
  if(dep[a] < dep[b])std::swap(a,b);
  int d = dep[a] - dep[b];
  for(int i = 17; i >= 0; i--) if(d>>i&1) a = f[i][a];
  if(a!=b) {
    for(int i = 17; i >= 0; i--) if(f[i][a]!=f[i][b]) a=f[i][a],b=f[i][b];
    a = f[0][a];
  }
  return a;
}
void lca_init() {
  for(int i = 0; i < 18; i++)for(int j = 1; j <= nn; j++) f[i][j] = 0;
  memset(dep,0,sizeof(dep));
  dep[0] = -1;
  memset(vis,false,sizeof(vis));
  for(int i = 1; i <= nn; i++) if(!vis[i]) DFS(i,0);
}
int main() {
  int u,v,q,m;
  int id = 0;
  int ca= 0;
  while(scanf("%d%d",&n,&m)!=EOF){
    ++ca;
    if(n==0&&m==0) break;
    E = 0;
    memset(head,-1,sizeof(head));
    for(int i = 0; i < m; i++){
      scanf("%d%d",&u,&v);
      add(u,v,i); add(v,u,i);
    }
    SCC();
    lca_init();
    scanf("%d",&q);
    while(q--){
      int id1,id2;
      scanf("%d%d",&id1,&id2);
      int a = belong[id1], b = belong[id2];
      printf("%d\n",(dep[a]+dep[b]-2*dep[LCA(a,b)])/2);
    }
  }
  return 0;
}

