#include <cstdio>
#include <ctime>
#include <vector>
#include <algorithm>
const int N = 100010;
const int M = N << 1;
int first[N], nxt[M], pnt[M], E;
void add_edge(int a, int b)
{
        pnt[E] = b;
        nxt[E] = first[a];
        first[a] = E++;
}
int parent[N], success[N], depth[N], size[N];
void dfs(int u, int fa=-1)
{
        parent[u] = fa, size[u] = 1, success[u] = -1;
        depth[u] = ~fa ? depth[fa] + 1 : 0;
        int v, mx = 0;
        for(int i = first[u]; ~i; i = nxt[i]) {
                if(fa == (v = pnt[i])) {
                        continue;
                }
                dfs(v, u);
                if(size[v] > mx) {
                        success[u] = v, mx = size[v];
                }
                size[u] += size[v];
        }
}
int head[N], id[N], cnt, mp[N];
void make_chain(int u, bool in_chain=false)
{
        head[u] = in_chain ? head[parent[u]] : u;
        id[u] = ++cnt;
        mp[cnt] = u;
        if(~success[u]) {
                make_chain(success[u], true);
        }
        for(int i = first[u]; ~i; i = nxt[i]) {
                if(success[u] == pnt[i] || parent[u] == pnt[i]) {
                        continue;
                }
                make_chain(pnt[i]);
        }
}
std::vector<int> operations[N];
void color(int a, int b, int c)
{
        while(head[a] != head[b]) {
                if(depth[head[a]] < depth[head[b]]) {
                        std::swap(a, b);
                }
                operations[id[head[a]]].push_back(c);
                operations[id[a] + 1].push_back(-c);
                a = parent[head[a]];
        }
        if(depth[a] > depth[b]) {
                std::swap(a, b);
        }
        operations[id[a]].push_back(c);
        operations[id[b] + 1].push_back(-c);
}
int ans[N];
struct Tree* pool;
struct Tree
{
        Tree *lc;
        Tree *rc;
        int mx, l, r;
        Tree(int l=0, int r=0): l(l), r(r) {
                mx = 0;
                if(l == r) {
                        return ;
                }
                int m = l + r >> 1;
                lc = new (pool++)Tree(l, m);
                rc = new (pool++)Tree(m + 1, r);
        }
        void insert(int p, int v) {
                if(l == r) {
                        mx += v;
                        return ;
                }
                int m = l + r >> 1;
                if(p <= m) {
                        lc->insert(p, v);
                } else {
                        rc->insert(p, v);
                }
                mx = std::max(lc->mx, rc->mx);
        }
        int query() {
                if(l == r) {
                        return l;
                }
                if(lc->mx == mx) {
                        return lc->query();
                } else {
                        return rc->query();
                }
        }
}node[1 << 18], *tree;
int max_color;
void solve(int n)
{
        pool = node;
        tree = new(pool++)Tree(1, max_color);
        for(int i = 1; i <= n; i++) {
                int sz = operations[i].size();
                for(int j = 0; j < sz; j++) {
                        int val = operations[i][j];
                        if(val > 0) {
                                tree->insert(val, 1);
                        } else {
                                tree->insert(-val, -1);
                        }
                }
                if(tree->mx != 0) {
                        ans[mp[i]] = tree->query();
                } else {
                        ans[mp[i]] = 0;
                }
        }
        for(int i = 1; i <= n; i++) {
                printf("%d\n", ans[i]);
        }
}
void init(int n)
{
        E = cnt = 0;
        std::fill(first, first + n + 1, -1);
        for(int i = 1; i <= n; i++) {
                operations[i].clear();
        }
}
bool read (int &x) {
        int c = getchar (); int sign = 1;
        while (~c && c < '0' || c > '9') { if (c == '-') sign = -1; c = getchar (); }
        for (x = 0; ~c && '0' <= c && c <= '9'; c = getchar ()) x = x * 10 + c - '0';
        x *= sign;
        return ~c;
}
int main()
{
        int n, q, a, b, c;
        while(scanf("%d%d", &n, &q) == 2) {
                if(!n && !q) break;
                init(n);
                for(int i = 2; i <= n; i++) {
                        a = i;
                        read(a), read(b);
                        add_edge(a, b), add_edge(b, a);
                }
                dfs(1), make_chain(1), parent[1]=0;
                max_color = 0;
                for(int i = 0; i < q; i++) {
                        read(a), read(b), read(c);
                        if(c > max_color)
                                max_color = c;
                        color(a, b, c);
                }
                if(q) {
                        solve(n);
                } else {
                        for(int i = 1; i <= n; i++) {
                                printf("0\n");
                        }
                }

        }
        return 0;
}
