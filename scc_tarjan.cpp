#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
const int N = 150010;
const int M = 300010;
struct Edge
{
        int head[N], nxt[M], pnt[M], E;
        void init() {
                E = 0;
                memset(head, -1, sizeof(head));
        }
        void add_edge(int a, int b){
                pnt[E] = b;
                nxt[E] = head[a];
                head[a] = E++;
        }
}ori, ne, rne;
int in[N], bel[N], dfn[N], low[N], stk[N], n, m, Time, top, btype;
std::vector<int> con[N];
void dfs(int u)
{
        low[u] = dfn[u] = ++Time; in[u] = 1; stk[++top] = u;
        for(int i = ori.head[u]; ~i; i = ori.nxt[i]) {
                int v = ori.pnt[i];
                if(!dfn[v]) {
                        dfs(v);
                        low[u]=std::min(low[u], low[v]);
                } else if(in[v]) {
                        low[u]=std::min(low[u], low[v]);
                }
        }
        if(low[u] == dfn[u]) {
                btype++; int v;
                do {
                        v = stk[top--], in[v] = 0, bel[v] = btype;
                        con[btype].push_back(v);
                }while(v != u);
        }
}
void scc()
{
        std::fill(in, in + n + 1, 0);
        std::fill(dfn, dfn + n + 1, 0);
        btype = Time = top = 0;
        for(int i = 1; i <= n; i++) if(!dfn[i]) {
                dfs(i);
        }
}
