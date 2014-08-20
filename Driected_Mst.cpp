#include <cstdio>
#include <algorithm>
typedef int item;
const int inf = ~0u>>2;
const int N = 1010;
const int M = 3010;
struct Edge
{
        int u, v;
        int cost;
}edge[M];
int pre[N], id[M], f[N], e;
item in[N];
void add_edge(int a, int b, int cost)
{
        edge[e].u = a; edge[e].v = b; edge[e].cost = cost;
        e++;
}
item directed_mst(int root, int n, int m)
{
        item ret = 0;
        int u, v;
        while(1) {
                std::fill(in, in + n, inf);
                for(int i = 0; i < m; i++) {
                        u = edge[i].u, v = edge[i].v;
                        if(edge[i].cost < in[v] && u != v) {
                                pre[v] = u;
                                in[v] = edge[i].cost;
                        }
                }
                for(int i = 0; i < n; i++) {
                        if(i != root && in[i] == inf) {
                                return -1;
                        }
                }
                int cirs = 0;
                std::fill(id, id + n, -1);
                std::fill(f, f + n, -1);
                in[root] = 0;
                for(int i = 0; i < n; i++) {
                        ret += in[i];
                        int v = i;
                        while(f[v] != i && id[v] == -1 && v != root) {
                                f[v] = i;
                                v = pre[v];
                        }
                        if(v != root && id[v] == -1) {
                                for(int u = pre[v]; u != v; u = pre[u]) {
                                        id[u] = cirs;
                                }
                                id[v] = cirs++;
                        }
                }
                if(cirs == 0) {
                        break;
                }
                for(int i = 0; i < n; i++) {
                        if(id[i] == -1) {
                                id[i] = cirs++;
                        }
                }
                for(int i = 0; i < m; i++) {
                        v = edge[i].v;
                        edge[i].u = id[edge[i].u];
                        edge[i].v = id[edge[i].v];
                        if(edge[i].u != edge[i].v) {
                                edge[i].cost -= in[v];
                        }
                }
                n = cirs;
                root = id[root];
        }
        return ret;
}
