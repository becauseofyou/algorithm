#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;
typedef int item;

#define rep(i,n) for(int i=0;i<n;i++)
const int inf = ~0u>>2;
const int M = 100010;
const int N = 100010;
struct EDGE {
   int u,v;
   item cost;
}E[M];
int pre[N] , ID[M] , f[N];
item in[N];
item Directed_Mst(int root,int NV,int NE) {
    item ret = 0; int u , v;
    while(true) {
        fill(in,in+NV,inf);
        rep(i,NE) {
            u = E[i].u; v = E[i].v;
            if(E[i].cost < in[v] && u != v) {
                pre[v] = u;
                in[v] = E[i].cost;
            }
        }	
        rep(i,NV) {
            if(i == root) continue;
            if(in[i] == inf) return -1;
        }
        int circle = 0;
        rep(i,NV) ID[i] = f[i] = -1;
        in[root] = 0;
        rep(i,NV) {
            ret += in[i];
            int v = i;
            while(f[v] != i && ID[v] == -1 && v != root) {
                f[v] = i;
                v = pre[v];
            }
            if(v != root && ID[v] == -1) {
                for(int u = pre[v]; u != v; u = pre[u]) ID[u] = circle;
                ID[v] = circle ++;
            }
        }
        if(circle == 0) break;
        rep(i,NV) if(ID[i] == -1) ID[i] = circle ++;
        rep(i,NE) {
            v = E[i].v;
            E[i].u = ID[E[i].u];
            E[i].v = ID[E[i].v];
            if(E[i].u != E[i].v) {
                E[i].cost -= in[v];
            }
        }
        NV = circle;
        root = ID[root];
    }
	return ret;
}
