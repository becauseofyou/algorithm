/*
 * 结构体或类名：Max_Flow
 * 变量名: maxFlow;
 * 函数名 : max_flow()
 *
 */
const int N = 1100010;
const int M = 6100000;
const int INF = ~0u >> 2;
template<class T>
struct Max_Flow  {
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
        *tail++ = s; level[s] = 0;
        while(front < tail && level[t] == -1) {
            int u = *front++;
            for(int e = head[u]; e != -1; e = nxt[e]) {
                if(cap[e] > 0 && level[pnt[e]] < 0) {
                    level[pnt[e]] = level[u] + 1;
                    *tail ++ = pnt[e];
                }
            }
        }
        return level[t] != -1;
    }
    void push(int t, int &flow) {
        int mi = INF, p = pre[t];
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
    void dfs(int u, int t, int &flow) {
        if(u == t) {
            push(t, flow);
            return ;
        }
        for(int &e = cur[u]; e != -1; e = nxt[e]) {
            if(cap[e] > 0 && level[u] < level[pnt[e]]) {
                pre[pnt[e]] = e;
                dfs(pnt[e], t, flow);
                if(level[sign] < level[u]) {
                    return ;
                }
                sign = t;
            }
        }
    }
    int dinic() {
        pre[s] = -1;
        int flow = 0;
        while(bfs(s, t)) {
            copy(head, head + n, cur);
            dfs(s, t, flow);
        }
        return flow;
    }
    int s, t, n;
    int Q[N], sign;
};
