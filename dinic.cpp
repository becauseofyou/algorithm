/*
 * 结构体或类名：Max_Flow
 * 变量名: maxFlow;
 * 函数名 : max_flow()
 *
 */
#include <cstdio>
#include <algorithm>
using std::fill;
using std::copy;
using std::min;
const int N = 20010;
const int M = 500010;
const int INF = ~0u >> 2;
template<class T>
struct Max_Flow {
    int s, t, n;
    int Q[N], sign;
    int head[N], level[N], cur[N], pre[N];
    int nxt[M], pnt[M], E;
    T cap[M];
    void init(int n, int s, int t) {
        this->n = n;
        this->s = s;
        this->t = t;
        E = 0;
        fill(head, head + n, -1);
    }
    void add_edge(int from, int to, T c) {
        pnt[E] = to;
        cap[E] = c;
        nxt[E] = head[from];
        head[from] = E++;
    }
    bool bfs(int s, int t) {
        sign = t;
        fill(level, level + n, -1);
        int *front = Q, *tail = Q;
        *tail++ = t; level[t] = 0;
        while(front < tail && level[s] == -1) {
            int u = *front++;
            for(int e = head[u]; e != -1; e = nxt[e]) {
                if(cap[e ^ 1] > 0 && level[pnt[e]] < 0) {
                    level[pnt[e]] = level[u] + 1;
                    *tail ++ = pnt[e];
                }
            }
        }
        return level[s] != -1;
    }
    void push(int t, T &flow) {
        T mi = INF;
        int p = pre[t];
        for(int p = pre[t]; p != -1; p = pre[pnt[p ^ 1]]) {
            mi = min(mi, cap[p]);
        }
        for(int p = pre[t]; p != -1; p = pre[pnt[p ^ 1]]) {
            cap[p] -= mi;
            if(!cap[p]) {
                sign = pnt[p ^ 1];
            }
            cap[p ^ 1] += mi;
        }
        flow += mi;
    }
    void dfs(int u, int t, T &flow) {
        if(u == t) {
            push(t, flow);
            return ;
        }
        for(int &e = cur[u]; e != -1; e = nxt[e]) {
            if(cap[e] > 0 && level[u] - 1 == level[pnt[e]]) {
                pre[pnt[e]] = e;
                dfs(pnt[e], t, flow);
                if(level[sign] > level[u]) {
                    return ;
                }
                sign = t;
            }
        }
    }
    T dinic() {
        pre[s] = -1;
        T flow = 0;
        while(bfs(s, t)) {
            copy(head, head + n, cur);
            dfs(s, t, flow);
        }
        return flow;
    }
};

int main() {
    int m, a, b, w, n;
    while(scanf("%d%d",&n,&m) != EOF) {
        Max_Flow<int> *gao = new Max_Flow<int>();
        gao->init(n + 2, 0, n + 1);
        int s = 0, t = n + 1;
        for(int i = 1; i <= n; i++) {
            scanf("%d%d", &a, &b);
            gao->add_edge(s, i, a);
            gao->add_edge(i, s, a);
            gao->add_edge(i, t, b);
            gao->add_edge(t, i, b);
        }
        for(int i = 0; i < m; i++) {
            scanf("%d%d%d",&a,&b,&w);
            gao->add_edge(a, b, w);
            gao->add_edge(b, a, w);
        }
        printf("%d\n",gao->dinic());
    }
    return 0;
}
