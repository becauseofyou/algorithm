#include <vector>
#include <cstdio>
#include <algorithm>
#include <cstring>
#define foreach(it,G) for(__typeof(G.begin())it = G.begin(); it != G.end(); it++)
const int N = 100010;
std::vector<int> G[N], rG[N], vs;
bool used[N];
int comp[N];
void add_edge(int a,int b) {
    G[a].push_back(b);
    rG[b].push_back(a);
}
void dfs(int u) {
    used[u] = true;
    foreach(it, G[u]) if(!used[*it]) dfs(*it);
    vs.push_back(u);
}
void rdfs(int u,int which) {
    used[u] = true;
    comp[u] = which;
    foreach(it,rG[u]) if(!used[*it]) rdfs(*it, which);
}
int scc(int n) {
    std::fill(used, used + n, false);
    vs.clear();
    for(int v = 0; v < n; v++) if(!used[v]) dfs(v);
    std::fill(used, used + n, false);
    int cnt = 0;
    for(int i = vs.size() - 1; i >= 0; i--)
        if(!used[vs[i]]) rdfs(vs[i],cnt++);
    return cnt;
}
